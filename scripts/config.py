def proteins():
    return ["PR", "IN", "RT"]


def algorithms():
    #return ["AlphaFold2", "AlphaFold3", "ESMFold", "ESM3-Open", "ESM3-Large", "Ember3D"]
    return ["AlphaFold2", "AlphaFold3", "ESMFold"]
# Note: Other algorithms (ESM3-Open, ESM3-Large, Ember3D) are currently excluded
# because they appear to have incomplete coverage for Protease (PR) compared to the PDB. --> I can't find a source for this.
