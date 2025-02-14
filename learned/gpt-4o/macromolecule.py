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
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt} < 1000 g/mol)"
    
    # Check for repeating structural motifs using substructure search (simplified approximation)
    # Typical motifs for proteins or glycosylated compounds and synthetic polymers
    repeating_motif_patterns = [
        "[NX3][CX3](=O)C", # Amide bonds (peptide linkage)
        "C(=O)O",         # Ester bonds (polyester)
        "[OX2H]",         # Hydroxyls, capping extended chains like in polysaccharides
    ]
    
    matches_count = sum(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in repeating_motif_patterns)
    if matches_count >= 2:
        # Check for complexity: high number of atoms and bonds
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        if num_atoms > 100 and num_bonds > 150:
            return True, f"High complexity and repeating structural motifs detected with {num_atoms} atoms and {num_bonds} bonds"
        
        return True, "Repeating structural motifs indicating a likely macromolecule"
    
    return False, "Structure does not meet macromolecule criteria"