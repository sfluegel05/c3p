"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: An alkanesulfonate oxoanion is one in which the sulfonate group –S(=O)(=O)[O–]
is directly attached to an sp³ (non‐aromatic, acyclic) carbon. In the original definition this carbon 
was meant to be at position 1 (the “alkyl” terminus) but some valid examples (e.g. isopropyl sulfonate)
show that the carbon may have two alkyl substituents provided those substituents are not aromatic.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    For our purposes we require that the molecule contains a –S(=O)(=O)[O–] group that is directly 
    attached to a carbon that is:
       • sp³ hybridized,
       • not in a ring,
       • and while it may have one or two carbon neighbors, any such carbon neighbor must not be aromatic.
       
    For example, taurine and its acyl derivatives or even isopropyl sulfonate (propane-2-sulfonate) would be acceptable.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alkanesulfonate oxoanion, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a sulfonate moiety attached to a carbon.
    # We start with a general pattern: a non‐aromatic carbon ([#6;!a]) directly attached to a sulfur 
    # that is doubly bonded to two oxygens and singly bonded to an oxygen with negative charge.
    sulfonate_pattern = Chem.MolFromSmarts("[#6;!a]S(=O)(=O)[O-]")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate pattern (C-S(=O)(=O)[O-]) found"
    
    # Iterate over each matching candidate; the match is a tuple of atom indices corresponding to:
    # (candidate carbon, sulfur, oxygen, oxygen, oxygen).
    for match in matches:
        candidate_idx = match[0]
        sulfur_idx = match[1]
        candidate = mol.GetAtomWithIdx(candidate_idx)
        sulfur = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check 1: candidate carbon must be sp3 hybridized.
        if candidate.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        # Check 2: candidate carbon must not be in a ring.
        if candidate.IsInRing():
            continue
        
        # Check 3: Determine the number of carbon neighbors (ignoring the bonded sulfur).
        # In our revised criteria, we allow one or two carbon neighbors.
        carbon_neighbors = []
        for neighbor in candidate.GetNeighbors():
            if neighbor.GetIdx() == sulfur.GetIdx():
                continue  # skip the attached sulfur
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors.append(neighbor)
        
        # Allow candidate carbon with at most 2 carbon neighbors.
        if len(carbon_neighbors) > 2:
            continue
        
        # Check 4: Ensure that any directly attached carbon neighbor is not aromatic;
        # if one is aromatic, then the local environment is not a simple alkyl chain.
        skip = False
        for nbr in carbon_neighbors:
            if nbr.GetIsAromatic():
                skip = True
                break
        if skip:
            continue
        
        # If all criteria are satisfied for this matching candidate, we classify this molecule as valid.
        return True, "Contains an alkanesulfonate oxoanion moiety attached to an sp3 (acyclic, non-aromatic) carbon"
    
    # If no valid candidate passed all filters, then we do not classify this molecule as an alkanesulfonate oxoanion.
    return False, "Found sulfonate group but not attached to a suitable sp3 (acyclic, non-aromatic) carbon"

# (Optional) Main block to run a few tests:
if __name__ == "__main__":
    test_smiles = [
        "C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCC",      # N-icosanoyltaurine(1-) should pass.
        "CC(=O)NCCCS([O-])(=O)=O",                        # acamprosate(1-) should pass.
        "OCCS([O-])(=O)=O",                              # isethionate should pass.
        "OC[C@H](O)CS([O-])(=O)=O",                       # (2S)-3-sulfopropanediol(1-) should pass.
        "CC(C)S([O-])(=O)=O",                            # propane-2-sulfonate; our revised method accepts a secondary (non-terminal) alkyl.
        "[K+].[K+].[H]C(=CC([H])=CC([H])=C1N(CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)C=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)S([O-])(=O)=O",  # False positive dye expected to fail.
    ]
    for sm in test_smiles:
        res, reason = is_alkanesulfonate_oxoanion(sm)
        print(f"SMILES: {sm}\nResult: {res}, Reason: {reason}\n")