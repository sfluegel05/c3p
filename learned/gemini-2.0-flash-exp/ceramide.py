"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides contain a sphingosine or sphinganine backbone, an amide-linked fatty acid,
    and may have a hydroxyl on the second carbon of the fatty acid. 

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a ceramide, False otherwise
         str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Sphingosine/Sphinganine Backbone identification (C-C-C with -OH,-NH2, -CH2OH)

    # Define SMARTS pattern for sphinganine and sphingosine core
    # Includes both saturated and unsaturated forms
    sphingosine_pattern = Chem.MolFromSmarts("[CX4H2][CX4H]([OX2H1])[CX4H]([NX3H2])([CX4H2][OX2H1])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
            return False, "No sphingosine/sphinganine backbone found."


    # 2. Amide Linkage Identification (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found."

    # 3. Fatty Acid Chain (long chain attached to carbonyl of amide linkage)
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found."
        
    
    # Check for reasonably long fatty acid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
      return False, "Fatty acid chain is too short."


    return True, "Contains sphingosine/sphinganine backbone with amide-linked fatty acid chain"