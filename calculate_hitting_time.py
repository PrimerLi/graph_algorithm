import os
import numpy as np
from scipy import *
from scipy.sparse import *
from scipy.sparse.linalg import inv
import time
import datetime

def read_transition_matrix(transition_matrix_file_name):
    assert(os.path.exists(transition_matrix_file_name))
    row_indices = []
    col_indices = []
    data = []
    reader = open(transition_matrix_file_name)
    for (index, string) in enumerate(reader):
        a = string.strip("\n").split(",")
        row = int(a[0])
        col = int(a[1])
        val = float(a[2])
        row_indices.append(row)
        col_indices.append(col)
        data.append(val)
    reader.close()
    dimension = max(max(row_indices), max(col_indices)) + 1
    return csr_matrix((data, (row_indices, col_indices)), shape = (dimension, dimension))

def get_mean_and_variance(transition_matrix):
    start = time.time()
    row, col = transition_matrix.shape
    assert(row == col)
    dimension = row
    row_indices = range(dimension)
    col_indices = [0] * dimension
    data = [1] * dimension
    ones = csr_matrix((data, (row_indices, col_indices)), shape = (dimension, 1))
    eps = 1.0e-16
    expectation = ones
    variance = csr_matrix(([0] * dimension, (row_indices, col_indices)), shape = (dimension, 1))
    power = ones
    expansion_order = 1000000
    counter = 0
    while counter < expansion_order:
        counter += 1
        power = transition_matrix * power
        expectation += power
        variance += counter * power
        error = np.linalg.norm(counter * power.todense())
        if True and counter % 500 == 0:
            print "Counter = " + str(counter) + ", error = " + str(error)
        if error < eps:
            break
    print "Counter = " + str(counter) + ", error = " + str(error)
    variance = 2.0 * variance
    variance = variance + expectation - expectation.multiply(expectation)
    end = time.time()
    print "Total time used in get_mean_and_variance = " + str(end - start) + " seconds. "
    #return inv(csr_matrix(([1] * dimension, (range(dimension), range(dimension))), shape = (dimension, dimension)) - transition_matrix) * ones, variance
    return expectation, variance

def read_vertices(vertices_file_name):
    assert(os.path.exists(vertices_file_name))
    vertices = []
    reader = open(vertices_file_name, "r")
    for (index, string) in enumerate(reader):
        vertices.append(string.strip("\n"))
    reader.close()
    return vertices

def read_adherents(adherents_file_name):
    assert(os.path.exists(adherents_file_name))
    reader = open(adherents_file_name, "r")
    adherents = []
    for (index, string) in enumerate(reader):
        if index == 0:
            target = string.strip("\n").split(" = ")[1]
        else:
            adherents.append(string.strip("\n"))
    reader.close()
    return target, adherents

def main():
    import sys
    if (len(sys.argv) != 4):
        print "transition_matrix_file_name = sys.argv[1], vertices_file_name = sys.argv[2], adherents_file_name = sys.argv[3]. "
        return -1

    transition_matrix_file_name = sys.argv[1]
    vertices_file_name = sys.argv[2]
    adherents_file_name = sys.argv[3]
    assert(os.path.exists(transition_matrix_file_name))
    assert(os.path.exists(vertices_file_name))
    assert(os.path.exists(adherents_file_name))

    vertices = read_vertices(vertices_file_name)
    transition_matrix = read_transition_matrix(transition_matrix_file_name)
    expectation, variance = get_mean_and_variance(transition_matrix)
    expectation = np.asarray(expectation.todense())[:, 0]
    variance = np.asarray(variance.todense())[:, 0]

    target, adherents = read_adherents(adherents_file_name)
    print "Saving expectation and variance results ... "
    assert(len(expectation) == len(variance))
    lines = []
    distances = dict()
    for i in range(len(expectation)):
        if vertices[i] == target:
            continue
        if vertices[i] in adherents:
            lines.append(vertices[i] + " -> " + target + ":" + str(1.0) + "," + str(0.0))
            distances[vertices[i]] = 1.0
        else:
            lines.append(vertices[i] + " -> " + target + ":" + str(expectation[i]) + "," + str(variance[i]))
            distances[vertices[i]] = expectation[i]
    
    sorted_lines = sorted(lines, key = lambda line: float(line.split(":")[-1].split(",")[0]))
    sorted_distances = []
    writer = open("hitting_times_mean_variance.txt", "w")
    for line in sorted_lines:
        writer.write(line + "\n")
        distance = float(line.split(":")[-1].split(",")[0])
        if distance == 1:
            continue
        sorted_distances.append(distance)
    writer.close()
   
    writer = open("sorted_distances.txt", "w")
    slopes = []
    for i in range(len(sorted_distances)):
        writer.write(str(i) + "  " + str(sorted_distances[i]) + "\n")
        if i < len(sorted_distances) - 1:
            slopes.append(sorted_distances[i+1] - sorted_distances[i])
    writer.close()

    max_slope = slopes[0]
    max_index = 0
    for i in range(1, len(slopes)):
        if max_slope < slopes[i]:
            max_slope = slopes[i]
            max_index = i

    within_community_vertices = [] + adherents
    distance_transition_point = sorted_distances[max_index+1]
    all_vertices = distances.keys()
    for vertex in all_vertices:
        if distances[vertex] < distance_transition_point:
            within_community_vertices.append(vertex)

    writer = open("within_community_vertices.txt", "w")
    writer.write("target = " + target + "\n")
    writer.write(",".join(within_community_vertices) + "\n")
    writer.close()

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
