def proteins():
    return ["PR", "IN", "RT"]


def algorithms():
    return ["AlphaFold2", "AlphaFold3", "ESMFold", "ESM3-Open", "ESM3-Large", "Ember3D"]

#TODO: move to plot
# start, end, kind
STRUCT = {
    "PR": [(13,19,'sheet'), (25,28,'helix'), (45,55,'sheet'), (59,64,'sheet'),
           (67,72,'sheet'), (76,81,'sheet'), (84,90,'sheet'), (95,99,'helix')],
    "IN": [(11, 23, 'helix'), (32, 42, 'helix'), (55, 63, 'sheet'), (70, 82, 'helix'),
        (101, 112, 'helix'), (121, 133, 'sheet'), (138, 147, 'sheet')],
    "RT": [(14, 24, 'helix'), (38, 47, 'helix'), (57, 66, 'sheet'), (71, 81, 'sheet'),
        (94, 105, 'helix'), (119, 128, 'helix'), (145, 154, 'sheet'), (174, 183, 'helix'),
        (186, 191, 'sheet'), (194, 200, 'sheet'), (204, 209, 'helix'),
        (212, 220, 'sheet'), (229, 238, 'helix'), (241, 246, 'sheet'),
        (250, 259, 'helix'), (265, 269, 'sheet'), (273, 278, 'helix'),
        (283, 295, 'helix'), (298, 304, 'sheet'), (312, 321, 'helix'),
        (333, 339, 'helix'), (342, 349, 'sheet'), (351, 355, 'sheet'),
        (361, 368, 'helix'), (371, 378, 'helix'), (382, 390, 'helix'),
        (399, 408, 'helix'), (410, 416, 'sheet'), (423, 430, 'sheet'),
        (437, 448, 'helix'), (454, 461, 'sheet'), (468, 479, 'helix'),
        (487, 494, 'helix'), (500, 506, 'sheet'), (512, 519, 'helix'),
        (523, 529, 'helix'), (534, 543, 'sheet'), (547, 553, 'sheet')],
}
