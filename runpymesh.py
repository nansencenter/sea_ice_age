from argparse import ArgumentParser
from math import radians

import numpy as np
import pymesh
from PyMesh import ShortEdgeRemoval, IsolatedVertexRemoval, FinFaceRemoval, LongEdgeRemoval, ObtuseTriangleRemoval

parser = ArgumentParser()
parser.add_argument('min_edge', type=float)
parser.add_argument('max_edge', type=float)
parser.add_argument('-f', '--file', type=str, default='mesh1.npz')
args = parser.parse_args()

filename = args.file
min_edge = args.min_edge
max_edge = args.max_edge

x = np.load(filename)['x']
y = np.load(filename)['y']
faces = np.load(filename)['t']
m = np.load(filename)['m']
z = np.zeros(x.size)
vertices = list(zip(x,y,z))
shapes_input = (x.size, faces.shape[0])
mesh = pymesh.form_mesh(vertices, faces)
important = np.union1d(mesh.boundary_vertices, m)
importance = np.zeros(mesh.vertices.shape[0], int)
importance[important] = -10

for i in range(5):
    importance_3  = importance[faces]
    remover = LongEdgeRemoval(vertices, faces)
    remover.run(max_edge)
    vertices = remover.get_vertices()
    faces = remover.get_faces()
    ori_faces = remover.get_ori_faces()
    importance_3 = importance_3[ori_faces]
    important = np.unique(faces[importance_3 == -10])
    importance = np.zeros(vertices.shape[0], int)
    importance[important] = -10

    collapser = ShortEdgeRemoval(vertices, faces)
    collapser.set_importance(importance)
    num = collapser.run(min_edge)
    vertices = collapser.get_vertices()
    faces = collapser.get_faces()
    face_indices = collapser.get_face_indices()
    importance_3 = importance_3[face_indices]
    important = np.unique(faces[importance_3 == -10])
    importance = np.zeros(vertices.shape[0], int)
    importance[important] = -10

    remover = FinFaceRemoval(vertices, faces)
    num = remover.run()
    vertices = remover.get_vertices()
    faces = remover.get_faces()
    face_indices = remover.get_face_indices()
    importance_3 = importance_3[face_indices]
    important = np.unique(faces[importance_3 == -10])
    importance = np.zeros(vertices.shape[0], int)
    importance[important] = -10

    remover = IsolatedVertexRemoval(vertices, faces)
    num = remover.run()
    vertices = remover.get_vertices()
    faces = remover.get_faces()
    ori_vertex_indices = remover.get_ori_vertex_indices()
    importance = importance[ori_vertex_indices]

xo,yo,zo = vertices.T
mo = np.where(importance == -10)[0]
np.savez(f'mesh2.npz', x=xo, y=yo, t=faces, m=mo)
print('Mesh clean up:', shapes_input, (xo.size, faces.shape[0]))
