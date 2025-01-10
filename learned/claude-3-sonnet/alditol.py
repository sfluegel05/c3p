"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:15972 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alditol(smiles: str):
    """
    Determines if a molecule contains an alditol moiety based on its SMILES string.
    An alditol is formally derivable from an aldose by reduction of the carbonyl group,
    having the general pattern HOCH2[CH(OH)]nCH2OH, where the hydroxyls may be substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an alditol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First look for basic alditol backbone pattern:
    # A chain of carbons where each carbon has an oxygen attached
    # The SMARTS pattern allows for substituted oxygens
    alditol_backbone = """
        # Start with carbon-oxygen
        [CX4]([OX2])
        # Connected to at least 2 more similar carbons
        [CX4]([OX2])[CX4]([OX2])
        # End with carbon-oxygen
        [CX4]([OX2])
    """
    alditol_pattern = Chem.MolFromSmarts(alditol_backbone.replace('\n',''))
    
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "No alditol backbone found (needs continuous chain of carbons with oxygen substituents)"

    # Get the matches
    matches = mol.GetSubstructMatches(alditol_pattern)
    
    for match in matches:
        # Convert match atoms to a submolecule for analysis
        match_atoms = set(match)
        
        # Check if these carbons form a continuous chain
        # by looking at bonds between matched atoms
        is_continuous = True
        for i in range(len(match)-1):
            bond = mol.GetBondBetweenAtoms(match[i], match[i+1])
            if bond is None:
                is_continuous = False
                break
                
        if not is_continuous:
            continue
            
        # For each carbon in the match, verify it has proper connectivity
        valid_carbons = True
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Each carbon should:
            # - Be sp3 hybridized
            # - Have one oxygen substituent
            # - Have correct number of carbons attached (1 or 2 depending on position)
            if atom.GetHybridization() != Chem.HybridizationType.SP3:
                valid_carbons = False
                break
                
            # Count oxygen and carbon neighbors
            o_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)
            c_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
            
            if o_neighbors < 1:
                valid_carbons = False
                break
                
            # Terminal carbons should have 1 carbon neighbor
            # Internal carbons should have 2 carbon neighbors
            if atom_idx in (match[0], match[-1]):
                if c_neighbors != 1:
                    valid_carbons = False
                    break
            else:
                if c_neighbors != 2:
                    valid_carbons = False
                    break
        
        if valid_carbons:
            return True, "Contains valid alditol structure (continuous carbon chain with hydroxyl/substituted hydroxyl groups)"
            
    return False, "No valid alditol structure found"