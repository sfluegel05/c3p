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
    min_macromolecule_weight = 500  # Lower threshold for smaller macromolecules
    if mol_wt < min_macromolecule_weight:
        return False, f"Molecular weight {mol_wt:.2f} Da is less than the minimum macromolecule threshold {min_macromolecule_weight} Da"

    # Define more specific repeating unit patterns
    # Detect typical monomer links, e.g., peptide ($C(=O)N), glycosidic ($O-$C-O), or nucleic acid links
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@H]1C([C@@H](O)C(O)[C@@H]1O)>>O[C@@H]([C@H]2OC)[C@H](O)[C@H](O2)")
    
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    if len(peptide_matches) >= 3 or len(glycosidic_matches) >= 3:
        return True, "Molecule is considered a macromolecule based on repeating peptide or glycosidic bond structures."

    return False, "No sufficient repeating unit motifs commonly found in macromolecules detected."