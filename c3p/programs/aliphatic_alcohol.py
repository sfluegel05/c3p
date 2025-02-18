"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'
Improvement rationale:
 - We first disqualify molecules that contain free carboxylic acid/carboxylate groups.
 - For the remaining molecules, we examine each –OH group.
   The –OH must be attached to a saturated (sp³), non‐aromatic carbon.
 - Then, from that candidate carbon we perform a depth‐first search (DFS)
   on its neighboring carbons (allowed if they are not aromatic, not in a ring and not carbonyl)
   to measure the longest contiguous aliphatic chain.
 - If any such chain reaches at least 6 carbons in length the molecule is classified as an aliphatic alcohol.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    Here, an aliphatic alcohol is defined as having at least one hydroxyl (-OH) group attached to a saturated 
    (sp³), non‐aromatic carbon. Furthermore, that candidate carbon must be connected to a contiguous, open-chain 
    aliphatic (acyclic) region (which may include one or two isolated sp2 centers) of at least 6 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Preliminary check: disqualify compounds with free carboxylic acid or carboxylate groups.
    # We consider a carboxylic acid pattern: [CX3](=O)[O;H1] and also its deprotonated form: [CX3](=O)[O-]
    acid_smarts = [Chem.MolFromSmarts("[CX3](=O)[O;H1]"), 
                   Chem.MolFromSmarts("[CX3](=O)[O-]")]
    for smarts in acid_smarts:
        if mol.HasSubstructMatch(smarts):
            return False, "Molecule contains carboxylic acid/carboxylate functionality"
    
    # Helper: Identify if a carbon can be considered part of a contiguous aliphatic chain.
    def is_allowed_chain_carbon(atom):
        # Must be a carbon that is not aromatic and not in any ring.
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing():
            return False
        # Exclude carbons that are carbonyl-like: any double bond to an oxygen.
        for bond in atom.GetBonds():
            # bond.GetBondTypeAsDouble() equals 2.0 for a double bond.
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return False
        return True
    
    # DFS to determine the maximum length (number of carbons) in the contiguous chain.
    # We count each allowed carbon exactly once per path.
    def dfs_chain(atom, visited):
        visited.add(atom.GetIdx())
        max_length = 1  # count this atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and is_allowed_chain_carbon(nbr) and nbr.GetIdx() not in visited:
                branch_length = dfs_chain(nbr, visited.copy())
                max_length = max(max_length, 1 + branch_length)
        return max_length

    # Iterate over oxygen atoms that might be part of an -OH group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Check that this oxygen has at least one hydrogen (explicit or implicit)
        if atom.GetTotalNumHs() < 1:
            continue  # not a free hydroxyl
        
        # Check the neighbors of the O; we are looking for a carbon to which the -OH is attached.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # Require candidate carbon to be saturated (sp³) and non‐aromatic.
            # (We relax the “not in ring” requirement so that, for example, side‐chains off cyclic structures are examined.)
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or nbr.GetIsAromatic():
                continue
            
            # Now, from this candidate carbon, ignore the oxygen we came from and search among its other neighbors
            # for a contiguous aliphatic chain (acyclic portion) of at least 6 carbons.
            for subnbr in nbr.GetNeighbors():
                # Do not follow back the hydroxyl oxygen.
                if subnbr.GetIdx() == atom.GetIdx():
                    continue
                # We require that the chain (if present) is in an acyclic, non‐aromatic region.
                if subnbr.GetAtomicNum() == 6 and is_allowed_chain_carbon(subnbr):
                    chain_length = dfs_chain(subnbr, set())
                    if chain_length >= 6:
                        return True, ("Found -OH group attached to a sp3, non-aromatic carbon "
                                      "that is connected to an acyclic aliphatic chain (chain length = {})."
                                      .format(chain_length))
    return False, ("No qualifying -OH group found attached to a sp3, non-aromatic carbon with an adjoining "
                   "acyclic aliphatic chain of at least 6 carbons.")

# (Optional) Some testing examples.
if __name__ == "__main__":
    test_smiles = [
        # True positives (expected):
        "CCCCCCC(C)O",           # octan-2-ol
        "O=C(OC)/C=C/CC(O)CCCCCC(O)C",  # Cladosporester A
        "CCCCCC(O)CC(=O)CCc1ccc(O)c(O)c1",  # 1-(3,4-Dihydroxyphenyl)-5-hydroxy-3-decanone
        "OCCCCCC/C=C/CCCC",       # 7E-dodecen-1-ol
        "O=C(CCCCCCCCCCCCCCCCCCCCCO)CCCCCC",  # 28-Hydroxy-7-octacosanone
        # False positives (expected to be disqualified due to acid functionality):
        "CCCCCCCCCCCCCCCCCC(O)C([O-])=O",  # 2-hydroxyarachidate
    ]
    
    for sm in test_smiles:
        classification, reason = is_aliphatic_alcohol(sm)
        print("SMILES:", sm)
        print("Classification:", classification)
        print("Reason:", reason)
        print("-" * 60)