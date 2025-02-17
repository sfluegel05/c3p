"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'

An aldopentose can appear as an open‐chain sugar (with a free aldehyde group) or as a cyclic (hemiacetal) sugar.
For a true aldopentose we enforce:
  1. The molecule must have exactly 5 carbon atoms.
  2. The molecule must have exactly 5 oxygen atoms (consistent with C5H10O5).
  3. If a free aldehyde group (SMARTS "[CX3H1](=O)") is found, it is classified as an open‐chain aldopentose.
  4. Otherwise, if the molecule is cyclic (has at least one ring) and does not contain a lactone
     (ester) functionality (SMARTS "[CX3](=O)[O]"), it is considered a cyclized aldopentose.
"""

from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose is defined as a pentose sugar (C5H10O5) that has an aldehyde group at one end in
    its open-chain form, or is present as the cyclic (hemiacetal) form that is in equilibrium with the 
    open-chain form. To qualify:
      - The molecule must have exactly 5 carbon atoms.
      - The molecule must have exactly 5 oxygen atoms.
      - It must have an aldehyde group (SMARTS "[CX3H1](=O)") or be a ring structure without a lactone (ester) 
        pattern (SMARTS "[CX3](=O)[O]").
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is an aldopentose; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 5:
        return False, f"Number of carbon atoms is {carbon_count}; expected 5 for an aldopentose"
    
    # 2. Count oxygen atoms.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 5:
        return False, f"Number of oxygen atoms is {oxygen_count}; expected 5 for C5H10O5"
        
    # 3. Define SMARTS patterns for aldehyde and lactone.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    
    # 4. If the molecule has a free aldehyde, classify it immediately as open-chain aldopentose.
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Open-chain aldopentose: aldehyde group detected."
    else:
        # 5. For cyclic forms: Check if the molecule is cyclic.
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            # Avoid molecules with lactone (ester) patterns.
            if mol.HasSubstructMatch(lactone_pattern):
                return False, "Cyclic structure contains lactone functionality; not an aldopentose."
            return True, "Cyclized aldopentose: no free aldehyde but cyclic, in equilibrium with open-chain form."
        else:
            return False, "Molecule is acyclic and lacks an aldehyde group; does not meet aldopentose criteria."
    
# (Optional testing calls could be added here, but the code block is self-contained.)