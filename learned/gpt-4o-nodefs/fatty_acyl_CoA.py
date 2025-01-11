"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    Fatty acyl-CoA is characterized by having Coenzyme A moiety which includes a pattern with pantetheine unit and a nucleotide segment,
    as well as a fatty acyl chain typically attached via a thiol ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for key substructures indicating Coenzyme A
    # Search for common Coenzyme A substructure patterns
    coenzymeA_pattern = Chem.MolFromSmarts("NCCSC(=O)C")
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)")

    # Check for presence of Coenzyme A segments
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Coenzyme A moiety pattern not found"
    
    if not mol.HasSubstructMatch(nucleotide_pattern):
        return False, "Nucleotide pattern not found"

    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate pattern not found"

    # Check for fatty acyl part via thiol ester and long carbon chains
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Thiol ester pattern not found"

    # Check for carbon chain length to verify fatty acyl part
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms for a typical fatty acid chain"

    return True, "Contains structural features of a fatty acyl-CoA"