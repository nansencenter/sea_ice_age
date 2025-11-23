import os
from contextlib import contextmanager

from sqlalchemy import (
    create_engine,
    Column,
    String,
    Float,
    ForeignKey,
    Index,
)
from sqlalchemy.orm import declarative_base, relationship, sessionmaker

DB_PATH = os.path.join(os.path.dirname(__file__), "zarr_index.sqlite")

Base = declarative_base()

class File(Base):
    __tablename__ = "files"

    path = Column(String, primary_key=True)   # absolute path
    mtime = Column(Float, nullable=False)     # last modification time

    arrays = relationship("Array", back_populates="file", cascade="all, delete-orphan")

class Array(Base):
    __tablename__ = "arrays"

    path = Column(String, ForeignKey("files.path"), primary_key=True)
    array_name = Column(String, primary_key=True)

    file = relationship("File", back_populates="arrays")

# index to speed up "WHERE array_name = ?"
Index("idx_arrays_name", Array.array_name)
_engine = create_engine(f"sqlite:///{DB_PATH}", future=True, connect_args={"timeout": 30.0})
SessionLocal = sessionmaker(bind=_engine, expire_on_commit=False)

def init_db():
    Base.metadata.create_all(_engine)

@contextmanager
def get_session():
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()

def get_file_arrays(path: str):
    """
    Return list of arrays for a file path.
    """
    init_db()
    abs_path = os.path.abspath(path)
    with get_session() as session:
        db_file = session.get(File, abs_path)
        if db_file is None:
            return []
        else:
            return [array.array_name for array in db_file.arrays]
