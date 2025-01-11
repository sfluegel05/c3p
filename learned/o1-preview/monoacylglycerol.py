"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: CHEBI:35700 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which any one of the hydroxyl groups is esterified with an acyl group,
    and the remaining two hydroxyl groups can be either unmodified (OH) or alkyl ethers (OR).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbone pattern: [CH2]-[CH]-[CH2]
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Define ester bond pattern: ester oxygen connected to carbonyl carbon
    ester_bond_pattern = Chem.MolFromSmarts("[C](=O)[O][CH]")
    
    # Loop over each glycerol backbone match
    for match in glycerol_matches:
        ester_count = 0
        non_ester_oxygen_count = 0
        glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Check each carbon in the glycerol backbone
        for carbon in glycerol_atoms:
            has_oxygen = False  # Flag to check if carbon has oxygen attached
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        has_oxygen = True
                        # Check if oxygen is part of an ester linkage
                        is_ester = False
                        for oxygen_neighbor in neighbor.GetNeighbors():
                            if oxygen_neighbor.GetIdx() != carbon.GetIdx() and oxygen_neighbor.GetAtomicNum() == 6:
                                # Check if this carbon is a carbonyl carbon
                                for nbr_bond in oxygen_neighbor.GetBonds():
                                    if nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr_bond.GetOtherAtom(oxygen_neighbor).GetAtomicNum() == 8:
                                        is_ester = True
                                        break
                        if is_ester:
                            ester_count += 1
                        else:
                            non_ester_oxygen_count +=1
            if not has_oxygen:
                break  # This carbon does not have any oxygen attached
        
        # Check if current glycerol backbone matches monoacylglycerol criteria
        if ester_count == 1 and (non_ester_oxygen_count == 2 or non_ester_oxygen_count == 1):
            return True, "Contains glycerol backbone with one esterified hydroxyl group and two hydroxyl/ether groups"
    
    return False, "Does not meet criteria for monoacylglycerol"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:35700',
        'name': 'monoacylglycerol',
        'definition': 'A glyceride in which any one of the R groups (position not specified) is an acyl group while the remaining two R groups can be either H or alkyl groups.',
        'parents': ['CHEBI:25212', 'CHEBI:17855']
    },
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 140,
    'num_false_positives': 3,
    'num_true_negatives': 182500,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.9790209790209791,
    'recall': 0.9032258064516129,
    'f1': 0.9395973154362416,
    'accuracy': 0.999914885568877
}