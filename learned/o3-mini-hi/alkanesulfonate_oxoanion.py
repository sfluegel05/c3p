"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: An alkanesulfonate in which the carbon at position 1 (attached to the sulfonate group)
represents the terminal carbon of an alkyl chain. In other words, the molecule must contain a sulfonate
group –S(=O)(=O)[O-] directly attached to an sp3 (non‐aromatic, acyclic) carbon that has at most one 
other carbon neighbour.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    For our purposes we require that the molecule contains a –S(=O)(=O)[O-] group directly 
    attached to a carbon that is:
       • sp3 hybridized,
       • not in a ring,
       • and terminal in an alkyl chain (i.e. that carbon has at most one carbon neighbor).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alkanesulfonate oxoanion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for sulfonate group attached to a non-aromatic carbon.
    # "[C;!a]S(=O)(=O)[O-]" means:
    # • a carbon that is not aromatic,
    # • that is directly attached to a sulfur bearing two double-bonded oxygens and one negative oxygen.
    sulfonate_pattern = Chem.MolFromSmarts("[C;!a]S(=O)(=O)[O-]")
    
    # Find all matches for that pattern.
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate pattern (C-S(=O)(=O)[O-]) found"
    
    # Iterate over each match and apply extra filtering.
    # The match indices correspond to (carbon, sulfur, oxygen, oxygen, oxygen).
    for match in matches:
        # Get the candidate carbon (first index) and the attached sulfur (second index)
        carbon = mol.GetAtomWithIdx(match[0])
        sulfur = mol.GetAtomWithIdx(match[1])
        
        # Check 1: Carbon must be sp3 hybridized.
        if carbon.GetHybridization() != rdchem.HybridizationType.SP3:
            continue  # skip if not sp3
        
        # Check 2: Carbon must not be in a ring.
        if carbon.IsInRing():
            continue  # skip if in a ring
            
        # Check 3: Verify that aside from the sulfur bond, the carbon is terminal.
        # Count the number of carbon neighbors (atomic num 6) other than the sulfur.
        carbon_neighbors = 0
        for neighbor in carbon.GetNeighbors():
            # Ignore the sulfur directly attached.
            if neighbor.GetIdx() == sulfur.GetIdx():
                continue
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
                
        # For a terminal carbon, we expect at most one carbon neighbor.
        if carbon_neighbors > 1:
            continue  # not terminal on an alkyl chain
        
        # (Optional extra check: Make sure that none of the other heavy neighbors are aromatic.
        # This was part of the earlier filter; here we mainly focus on the terminal nature.)
        valid = True
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetIdx() == sulfur.GetIdx():
                continue
            if neighbor.GetAtomicNum() > 1 and neighbor.GetIsAromatic():
                valid = False
                break
        if not valid:
            continue
        
        # If we reached here, we have a valid alkanesulfonate oxoanion moiety.
        return True, "Contains an alkanesulfonate oxoanion moiety attached to a terminal sp3 carbon"
        
    # If no match passes all filters, then the moiety is not in an appropriate alkane context.
    return False, "Found sulfonate group but not attached to a terminal sp3 (acyclic, non-aromatic) carbon"
    
# (Optional) Main block to run simple tests:
if __name__ == "__main__":
    test_smiles = [
        "C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCC",      # N-icosanoyltaurine(1-): should be True
        "CC(=O)NCCCS([O-])(=O)=O",                        # acamprosate(1-): True
        "OCCS([O-])(=O)=O",                              # isethionate: True
        "[K+].[K+].[H]C(=CC([H])=CC([H])=C1N(CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)C=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)S([O-])(=O)=O",  # NIR-3 dye; likely False
    ]
    for sm in test_smiles:
        res, reason = is_alkanesulfonate_oxoanion(sm)
        print(f"SMILES: {sm}\nResult: {res}, Reason: {reason}\n")