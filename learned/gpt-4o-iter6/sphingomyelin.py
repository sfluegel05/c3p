"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_sphingomyelin(smiles: str):
    """Determines if a molecule is a sphingomyelin."""
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amide linkage with fatty acid (N-C=O)
    amide_pattern = Chem.MolFromSmarts("N[C@@]([C@H](O)COP([O-])(=O)OCC[N+](C)(C)C)=O")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Does not have the characteristic amide linkage with fatty acid"
    
    # Check for phosphorylcholine group
    phosphorylcholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "Does not have phosphorylcholine group"
    
    # Check for sphingoid base structure - long chain with double bond/s and hydroxyl groups
    # This is a simplification and may not cover all structural features of sphingoid bases
    base_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@@H]/C=C/CCCCCCCC")
    if not mol.HasSubstructMatch(base_pattern):
        return False, "Missing sphingoid base features"
    
    return True, "Molecule matches the core structural features of sphingomyelin"