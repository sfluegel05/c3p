"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Review and calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Lower threshold to allow smaller macro-units
        return False, f"Molecular weight possibly too low ({mol_wt} < 500 g/mol)"
    
    # Define repeating structural motifs typical in macromolecules
    repeating_motif_patterns = [
        "[NX3][CX3](=O)[O,N]", # Amide bonds (peptide linkages)
        "COC(=O)",            # Ester bonds (polyester components)
        "[OX2H][CX4]",        # Hydroxyl caps, potential for chains like polysaccharides
        "C1OC1",              # Ether linkage in sugars
        "[PO4]",              # Phosphate backbone in nucleic acids
    ]
    
    # Count matches of repeating motifs
    matches_count = sum(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in repeating_motif_patterns)
    if matches_count >= 3:
        # Check for complexity: higher number of atoms and potential repeating units
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        if num_atoms > 60 and num_bonds > 80:
            return True, f"Complex structure with repeating motifs detected: {num_atoms} atoms, {num_bonds} bonds"
        
        return True, "Repeating structural motifs detected, potential macromolecule structure"
    
    return False, "Structure does not meet macromolecule criteria"