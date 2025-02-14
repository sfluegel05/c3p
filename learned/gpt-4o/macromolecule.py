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
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:  # Set a higher threshold for macromolecules
        return False, f"Molecular weight too low ({mol_wt} < 1000 g/mol)"
    
    # Define repeating structural motifs typical in macromolecules
    repeating_motif_patterns = [
        "[NX3][CX3](=O)[O,N]", # Amide bonds (peptide linkages)
        "COC(=O)",            # Ester bonds (polyester components)
        "[OX2H][CX4]",        # Hydroxyl caps, potential for chains like polysaccharides
        "C1OC1",              # Ether linkage in sugars
        "[PO4]",              # Phosphate group
    ]
    
    # Count matches of repeating motifs
    matches_count = sum(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in repeating_motif_patterns)
    # Increase threshold for number of motifs
    if matches_count >= 5:
        # Check for complexity: higher number of atoms and bonds
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        if num_atoms > 100 and num_bonds > 150:
            return True, f"Complex structure with repeating motifs detected: {num_atoms} atoms, {num_bonds} bonds"
        
        return True, "Extensive repeating structural motifs detected, indicating potential macromolecule structure"
    
    return False, "Structure does not meet macromolecule criteria"