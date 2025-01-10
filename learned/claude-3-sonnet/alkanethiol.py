"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: A compound in which a sulfanyl group, -SH, is attached to an alkyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Must have at least one carbon and one sulfur
    if not all(any(atom.GetAtomicNum() == n for atom in mol.GetAtoms()) 
               for n in [6, 16]):  # C and S
        return False, "Missing required elements (C and S)"

    # Look for -SH (sulfanyl) group pattern more specifically
    sulfanyl_pattern = Chem.MolFromSmarts("[SH1]-[CX4]")  # Must be connected to sp3 carbon
    if not mol.HasSubstructMatch(sulfanyl_pattern):
        return False, "No sulfanyl (-SH) group attached to alkyl carbon"

    # Count sulfanyl groups
    sulfanyl_matches = len(mol.GetSubstructMatches(sulfanyl_pattern))
    if sulfanyl_matches > 3:  # Allow up to 3 thiol groups
        return False, "Too many sulfanyl groups"

    # Exclude aromatic compounds
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Exclude compounds with complex functional groups
    exclusion_patterns = [
        ("[CX3](=O)[OX2H1]", "Contains carboxylic acid"),
        ("[CX3](=O)[OX2][CX4]", "Contains ester group"),
        ("[NX3](~[!#1])(~[!#1])~[!#1]", "Contains tertiary amine"),
        ("[#7X2]=,:[#6]", "Contains imine or similar"),
        ("[#6X3]~[#6X3]~[#6X3]", "Contains conjugated system"),
        ("[S;X3,X4,X5,X6]", "Contains oxidized sulfur"),
        ("[BX3]", "Contains boron"),
        ("[SiX4]", "Contains silicon"),
        ("[S;R][S;R]", "Contains disulfide ring"),
        ("[C;R7,R8,R9,R10]", "Contains large rings")
    ]

    for pattern, reason in exclusion_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            return False, reason

    # Allow simple rings (up to 6-membered)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:
        return False, "Contains multiple rings"
    
    # Check molecular complexity
    if rdMolDescriptors.CalcNumRotatableBonds(mol) > 15:
        return False, "Too many rotatable bonds"
        
    if rdMolDescriptors.CalcNumAtomStereoCenters(mol) > 3:
        return False, "Too many chiral centers"

    # Get all sulfur atoms and verify their environment
    valid_sulfanyl = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            if atom.GetDegree() == 1:  # Terminal sulfur (as in -SH)
                carbon = atom.GetNeighbors()[0]
                if carbon.GetAtomicNum() == 6:  # Connected to carbon
                    if not carbon.GetIsAromatic():
                        valid_sulfanyl = True

    if not valid_sulfanyl:
        return False, "No valid sulfanyl group found"

    return True, "Contains sulfanyl group (-SH) attached to an alkyl group"