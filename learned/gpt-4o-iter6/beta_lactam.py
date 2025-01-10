"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    Beta-lactams have a four-membered ring which includes the amide nitrogen
    and the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern to match a four-membered beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1(=O)NCC1")  # Improved SMARTS pattern
    
    # Check if the molecule matches the beta-lactam pattern
    if mol.HasSubstructMatch(beta_lactam_pattern):
        # Further verify the identified ring is indeed four-membered
        ring_info = mol.GetRingInfo()
        for ring_atoms in ring_info.AtomRings():
            if len(ring_atoms) == 4:
                # Check if the ring atoms match the SMARTS query
                query = rdqueries.AtomRingTestQuery(4)
                if all(query.Match(mol.GetAtomWithIdx(atom)) for atom in ring_atoms):
                    return True, "Contains a four-membered beta-lactam ring"
    
    return False, "Does not contain a four-membered beta-lactam ring"