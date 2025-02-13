"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36737 quinone

A quinone is a compound having a fully conjugated cyclic dione structure, such as that of
benzoquinones, derived from aromatic compounds by conversion of an even number of -CH= groups
into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and
heterocyclic analogues are included).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generate tautomers
    tautomers = [mol]
    try:
        tautomers.extend(AllChem.EnumerateTautomers(mol))
    except ValueError:
        pass
    
    # Check each tautomer for quinone structure
    for tautomer in tautomers:
        # Look for conjugated cyclic dione pattern
        dione_patterns = [
            rdqueries.GetQueryForSubstructureMatch(Chem.MolFromSmarts("[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:1]=[O&R1]:[a&r]:[O&R1]:[a&r]1")),
            rdqueries.GetQueryForSubstructureMatch(Chem.MolFromSmarts("[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:[C&R1]=,:1]=[O&R1]:[a&r]:[O&R1]:[a&r]1"))
        ]
        if any(tautomer.HasSubstructMatch(pattern) for pattern in dione_patterns):
            # Check if carbonyl groups are part of an aromatic or conjugated system
            carboxyl_atoms = [atom for atom in tautomer.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0]
            if all(atom.IsInAromaticRing() or atom.IsInRingOfSize(5) or atom.IsInRingOfSize(6) for atom in carboxyl_atoms):
                return True, "Contains a fully conjugated cyclic dione structure derived from an aromatic compound"
    
    return False, "No quinone structure found"