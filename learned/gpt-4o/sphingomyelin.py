"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelin is characterized by a sphingoid base linked to a fatty acid
    via an amide bond and a phosphorylcholine headgroup.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more general SMARTS patterns for key structural features of sphingomyelin
    # General pattern for amide bond with any alkyl chain
    amide_bond_pattern = Chem.MolFromSmarts("[NX3][CX4]([H])[CX3](=O)[CX4]")
    
    # Phosphorylcholine group - verify if exactly the same structure needed
    phosphorylcholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    
    # General pattern for sphingoid base (ensuring variability in stereochemistry and functional groups)
    sphingoid_base_pattern = Chem.MolFromSmarts("C[C@@H](OC)C[CX4](O)(C)")

    # Check for amide bond
    if not mol.HasSubstructMatch(amide_bond_pattern):
        return False, "No amide bond with a fatty acid found"

    # Check for phosphorylcholine
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "No phosphorylcholine headgroup found"
    
    # Check for sphingoid base
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    return True, "Contains sphingoid base, amide-linked fatty acid, and phosphorylcholine headgroup"