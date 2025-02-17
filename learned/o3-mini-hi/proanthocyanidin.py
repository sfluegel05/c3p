"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.
Heuristic: We require that (a) the molecule is valid, (b) has a molecular weight above 500 Da,
(c) has at least 4 aromatic six‐membered rings (which would be the benzene rings from the flavan units),
and (d) has a high number of –OH groups (here at least 8). These filters are not perfect but can help as an initial classifier.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string using heuristic filters.
    We require:
      - a valid structure,
      - a molecular weight of at least 500 Da (to exclude monomers),
      - at least 4 aromatic rings of size 6 (the phenolic rings from two or more flavan units),
      - at least 8 hydroxyl groups (–OH) indicating polyhydroxylation.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a proanthocyanidin, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); proanthocyanidins are oligomers (expected >500 Da)"
    
    # Count aromatic benzene rings.
    # We get the rings from the molecule and then check for aromatic rings of size 6.
    ring_info = mol.GetRingInfo().AtomRings()
    benzene_ring_count = 0
    for ring in ring_info:
        if len(ring) == 6:
            # Check that every atom in the ring is aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                benzene_ring_count += 1
    if benzene_ring_count < 4:
        return False, f"Not enough aromatic benzene rings ({benzene_ring_count} found); expect at least 4 to suggest multi-unit flavan structure"
    
    # Count hydroxyl groups using a SMARTS for -OH.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_OH = len(oh_matches)
    if n_OH < 8:
        return False, f"Not enough hydroxyl groups ({n_OH} found); proanthocyanidins are typically highly hydroxylated"
    
    # If all criteria are met, we assume the structure is a proanthocyanidin.
    return True, ("Molecule has high molecular weight, multiple aromatic (benzene) rings, and many hydroxyl groups "
                  "consistent with being a proanthocyanidin (a flavonoid oligomer)")

# Note: This heuristic classifier may miss unusual cases or borderline molecules.