"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:38159 hopanoid

A hopanoid is a triterpenoid based on a hopane skeleton.
The hopane skeleton consists of a pentacyclic ring system
with 5 rings arranged in a particular way.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for hopane skeleton pattern
    hopane_pattern = Chem.MolFromSmarts("[C@@]12[C@]([C@@]3([C@]([C@@]4([C@]([C@@]5([C@H](C(CCC5)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1")
    if not mol.HasSubstructMatch(hopane_pattern):
        return False, "No hopane skeleton found"
    
    # Check number of rings and carbon skeleton
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 5:
        return False, "Too few rings for hopanoid"
    
    # Check molecular formula - hopanoids typically C30-C35
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    c_count = mol_formula.count("C")
    if c_count < 30 or c_count > 35:
        return False, f"Carbon count {c_count} outside typical range for hopanoids"
    
    # Additional checks based on chemical knowledge
    # ...
    
    return True, "Contains hopane skeleton and meets other criteria for hopanoids"

__metadata__ = {"chemical_class": {"id":"CHEBI:38159", "name":"hopanoid", "definition":"A triterpenoid based on a  hopane skeleton.","parents":["CHEBI:35708"]},
"config":{"llm_model_name":"lbl/claude-sonnet","f1_threshold":0.8,"max_attempts":5,"max_positive_instances":None,"max_positive_to_test":None,"max_negative_to_test":None,"max_positive_in_prompt":50,"max_negative_in_prompt":20,"max_instances_in_prompt":100,"test_proportion":0.1},
"message":"",
"attempt":0,
"success":True,
"best":True,
"error":"",
"stdout":None,
"num_true_positives":207,
"num_false_positives":0,
"num_true_negatives":182397,
"num_false_negatives":1,
"num_negatives":None,
"precision":1.0,
"recall":0.9951928963919886,
"f1":0.9975865747642091,
"accuracy":0.9999994560358548}