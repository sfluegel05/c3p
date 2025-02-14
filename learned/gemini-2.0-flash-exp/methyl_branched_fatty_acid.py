"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid has a carboxylic acid group, a hydrocarbon chain,
    and only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ring structures - not allowed
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[R]")):
        return False, "Molecule contains rings"

    # Check for non-methyl branches
    non_methyl_branch = Chem.MolFromSmarts("C[CH2,CH]([CH3])")
    if mol.HasSubstructMatch(non_methyl_branch):
         return False, "Molecule contains non-methyl branches"

    # Find main chain carbons connected to the acid group
    chain_carbon = mol.GetSubstructMatches(Chem.MolFromSmarts("CC(=O)O"))
    if len(chain_carbon) == 0:
        return False, "Chain not connected to acid"

    # Get number of carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)


    #check main chain
    
    main_chain_carbons = set()
    
    for match in chain_carbon:
      main_chain_start = match[0]
      main_chain_carbons.add(main_chain_start)

      
    #Find adjacent chain
    next_carbon_atoms = []
    for atom_id in main_chain_carbons:
        atom = mol.GetAtomWithIdx(atom_id)
        for neighbor in atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != chain_carbon[0][1]:
            next_carbon_atoms.append(neighbor.GetIdx())

    main_chain_carbons.update(next_carbon_atoms)



    # Check for minimum chain length of at least 4 carbons excluding the carbonyl carbon (C=O).
    if num_carbons < 4:
      return False, "Too few carbons to be a fatty acid"

    # Get number of branches
    branch_pattern = Chem.MolFromSmarts("[CH]([CH3])")
    num_branches = len(mol.GetSubstructMatches(branch_pattern))


    #Check if we have other functional groups
    other_functional_group_pattern = Chem.MolFromSmarts("[!#1!#6!#8!#7!#15]")
    if mol.HasSubstructMatch(other_functional_group_pattern):
        return False, "Molecule contains other functional groups"

    return True, "Methyl-branched fatty acid"