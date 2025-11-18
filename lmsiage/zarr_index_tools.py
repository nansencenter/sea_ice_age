import os
import glob

from .mesh_file import MeshFile
from .zarr_index import init_db, get_session, File, Array


def list_arrays_in_file(filepath: str):
    """Use MeshFile to list array names in one zip Zarr file."""
    mf = MeshFile(filepath)
    return mf.read_names()

def update_index_for_dir(root_dir: str, pattern: str = "*.zip"):
    """
    Scan root_dir recursively for zip Zarr files and index them.
    Re-indexes a file only if mtime changed.
    """
    init_db()
    root_dir = os.path.abspath(root_dir)

    with get_session() as session:
        for zip_path in glob.iglob(os.path.join(root_dir, "**", pattern), recursive=True):
            zip_path = os.path.abspath(zip_path)

            try:
                st = os.stat(zip_path)
            except FileNotFoundError:
                continue
            mtime = st.st_mtime

            # Check if file already indexed and unchanged
            db_file = session.get(File, zip_path)
            if db_file is not None and abs(db_file.mtime - mtime) < 1e-6:
                continue  # up-to-date

            # Read arrays from the file
            try:
                array_names = list_arrays_in_file(zip_path)
            except Exception as e:
                print(f"Failed to read {zip_path}: {e}")
                continue

            if db_file is None:
                db_file = File(path=zip_path, mtime=mtime)
                session.add(db_file)
            else:
                db_file.mtime = mtime
                # existing Array entries will be cleared by reassigning db_file.arrays

            # Replace arrays for that file
            db_file.arrays.clear()
            for name in array_names:
                db_file.arrays.append(Array(path=zip_path, array_name=name))

def files_missing_array(array_name: str):
    """
    Return list of file paths that do NOT contain given array_name.
    """
    init_db()
    with get_session() as session:
        # subquery: all paths that DO have this array
        from sqlalchemy import select
        subq = select(Array.path).where(Array.array_name == array_name) #.subquery()
        q = select(File.path).where(~File.path.in_(subq))
        res = session.execute(q).scalars().all()
    return res

def files_with_array(array_name: str):
    """
    Return list of file paths that DO contain given array_name.
    """
    init_db()
    with get_session() as session:
        from sqlalchemy import select
        q = (
            select(File.path)
            .join(Array, File.path == Array.path)
            .where(Array.array_name == array_name)
        )
        res = session.execute(q).scalars().all()
    return sorted(res)

def cleanup_missing_files():
    """
    Remove File (and related Array) entries for files that no longer exist on disk.
    """
    init_db()
    with get_session() as session:
        from sqlalchemy import select
        paths = session.execute(select(File.path)).scalars().all()
        missing = [p for p in paths if not os.path.exists(p)]
        if not missing:
            return
        # delete in bulk; arrays are removed via cascade
        session.query(File).filter(File.path.in_(missing)).delete(synchronize_session=False)

def reset_index():
    """
    Remove all entries from all tables in the zarr index database.
    """
    init_db()
    with get_session() as session:
        # Order matters because of FK: delete children (Array) first, then parents (File)
        session.query(Array).delete(synchronize_session=False)
        session.query(File).delete(synchronize_session=False)

def delete_file_from_index(path: str):
    """
    Remove a specific File and its associated Arrays from the index.

    Parameters
    ----------
    path : str
        Path to the file (absolute or relative). It will be normalized to
        an absolute path before lookup.
    """
    init_db()
    abs_path = os.path.abspath(path)
    with get_session() as session:
        db_file = session.get(File, abs_path)
        if db_file is None:
            return
        session.delete(db_file)  # cascades to Array via relationship

