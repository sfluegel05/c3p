"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid with one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for at least one hydroxyl group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(alcohol_pattern):
         return False, "No hydroxy group found"

    # Check for a long aliphatic carbon chain connected to the carboxylic acid
    # This pattern ensures that we select chain carbons (X4), attached to the carboxylic group C
    # It also allows for double bonds along the chain.
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)O[CX4][CX3,CX4]([CX4])~[CX3,CX4]~[CX3,CX4]~[CX3,CX4]~[CX3,CX4]")
    matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if not matches:
        return False, "No long carbon chain attached to carboxylic group found"

    # Check if at least ONE hydroxy group is on the fatty acid chain
    hydroxy_on_chain = False
    for match in matches:
      # get the carboxyl carbon, and the carbons of the chain
      chain_atoms = [match[3],match[4],match[5], match[6], match[7]]
      for atom in chain_atoms:
        # check if the chain atoms have an OH attached to them
        for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalValence() == 2:
                hydroxy_on_chain = True
                break
        if hydroxy_on_chain:
            break # Found one, we can exit both loops

    if not hydroxy_on_chain:
        # Check if the hydroxyl is on the chain *or* attached to the carboxyl group carbon
        carboxyl_carbon = mol.GetSubstructMatches(acid_pattern)[0][0]
        for neighbor in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalValence() == 2:
                return False, "Hydroxy not on fatty acid chain."
            

        #If the hydroxyl group is not on the chain, it's not a hydroxy fatty acid
        
        
    # Check if chain is long enough using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
      return False, "Carbon chain is too short"

    # Check the number of carbons:
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons to be a fatty acid"
    
    # Check if the carbon chain is mostly aliphatic (avoiding aromatics, cycles)
    # This is a weak check and might be improved by other methods
    num_sp2_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() != Chem.HybridizationType.SP3)
    if num_sp2_carbons > 2:
        return False, "Too many sp2 carbons, not a fatty acid"

    return True, "Contains a carboxylic acid group, at least one hydroxy group on the fatty acid chain and a long carbon chain"