"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:24116 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a molecule of high relative molecular mass, the structure
    of which essentially comprises the multiple repetition of units derived,
    actually or conceptually, from molecules of low relative molecular mass.

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
    
    # Macromolecules typically have a molecular weight > 1000 Da
    if mol_wt < 1000:
        return False, f"Molecular weight {mol_wt:.2f} Da is too low for a macromolecule"
    
    # Check for repetitive units (fingerprint-based)
    fp = Chem.RDKFingerprint(mol, maxPath=7, fpSize=2048)
    n_bits_set = len(fp.GetNonzeroElements())
    
    # Macromolecules tend to have a smaller number of unique substructures
    if n_bits_set > 100:
        return False, "Too many unique substructures for a macromolecule"
    
    # Look for polymeric repeat units (e.g., peptides, polysaccharides)
    repeat_units = ['C(=O)N', 'OCC', 'OC=O', 'NC=O']
    for pattern in repeat_units:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return True, "Contains repetitive units derived from smaller molecules"
    
    return None, None