"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analog based on its SMILES string.
    A nucleobase analog is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analog, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for pyrimidine or purine ring system
    pyrimidine_pattern = Chem.MolFromSmarts("c1[cn][cn][cn]1") # c1ncncn1
    purine_pattern = Chem.MolFromSmarts("c1[cn][cn][cn][cn]2c1[cn]2") # c1ncncn2c1cn2
    if not mol.HasSubstructMatch(pyrimidine_pattern) and not mol.HasSubstructMatch(purine_pattern):
        return False, "No pyrimidine or purine ring found."

    # 2. Check for common substituents
    # Define SMARTS patterns for common functional groups
    methyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    formyl_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]")
    amino_pattern = Chem.MolFromSmarts("[NX3H2]")
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    etheno_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")

    # Count occurrences
    num_methyl = len(mol.GetSubstructMatches(methyl_pattern))
    num_formyl = len(mol.GetSubstructMatches(formyl_pattern))
    num_hydroxy = len(mol.GetSubstructMatches(hydroxy_pattern))
    num_amino = len(mol.GetSubstructMatches(amino_pattern))
    num_carboxy = len(mol.GetSubstructMatches(carboxy_pattern))
    num_halogen = len(mol.GetSubstructMatches(halogen_pattern))
    num_carbonyl = len(mol.GetSubstructMatches(carbonyl_pattern))
    num_etheno = len(mol.GetSubstructMatches(etheno_pattern))

    # If at least one of them is found, then it is probably a nucleobase analog
    if num_methyl > 0 or num_formyl > 0 or num_hydroxy > 0 or num_amino > 0 or num_carboxy > 0 or num_halogen > 0 or num_carbonyl > 0 or num_etheno > 0:
       return True, "Contains a pyrimidine or purine ring with common substitutions."

    # 3. Handling Aza/Oxa analogs (we already check it is pyrimidine/purine)
    # Add here an extra check to identify nitrogen at specific positions
    azapurine_pattern = Chem.MolFromSmarts("n1[cn][cn][cn][cn]2c1[cn]2")
    azapyrimidine_pattern = Chem.MolFromSmarts("n1[cn][cn][cn]1")
    oxapurine_pattern = Chem.MolFromSmarts("o1[cn][cn][cn][cn]2c1[cn]2")
    oxapyrimidine_pattern = Chem.MolFromSmarts("o1[cn][cn][cn]1")


    if mol.HasSubstructMatch(azapurine_pattern) or mol.HasSubstructMatch(azapyrimidine_pattern) or mol.HasSubstructMatch(oxapurine_pattern) or mol.HasSubstructMatch(oxapyrimidine_pattern):
       return True, "Contains a substituted pyrimidine or purine ring."
    

    # If none of the above is satisfied, then its not a nucleobase analog
    return False, "Does not resemble a nucleobase."