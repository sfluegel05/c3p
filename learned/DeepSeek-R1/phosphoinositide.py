"""
Classifies: CHEBI:18179 phosphoinositide
"""
pi_core_pattern = Chem.MolFromSmarts(
    "[CH2]-[CH](-O-C(=O)-[!O])-[CH2]-O-P(=O)(-O-[C@]1-[C@H](O)-[C@H](O)-[C@H](O)-[C@H](O)-[C@H](O)-[C@H]1-O)-O"
)