"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
Revised classifier that relaxes the molecular weight threshold and improves core detection by requiring
that the key amino-diol (or carbonyl variant) motif is found on an acyclic fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds are defined as sphinganine (and its homologs/stereoisomers) and their
    hydroxy/unsaturated derivatives. These compounds typically have:
      - A long aliphatic (acyclic) chain (>=8 connected sp3 carbon atoms).
      - A characteristic sphingoid core featuring a 2-amino-1,3-diol or its carbonyl variant.
      - At least one nitrogen atom (from the amino group).
      - A molecular weight that usually is not extremely low (we set a lower bound near 200 Da).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as sphingoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for sufficient long alkyl chain.
    # We require a contiguous chain of at least 8 carbon atoms.
    # (This SMARTS looks for a chain "CCCCCCCC".)
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long aliphatic chain (>=8 carbons) found."
    
    # Check total carbon count. (Sphingoid compounds normally are long-chain molecules.)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms ({c_count}); sphingoid molecules typically have a long chain."
    
    # Require at least one nitrogen (for the amino group)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found; sphingoid compounds require an amino group."
    
    # Define SMARTS patterns for the typical sphingoid core motifs.
    # Pattern 1: A typical amino-diol motif: (ignoring stereochemistry) C(O)C(N)CO
    core_pattern1 = Chem.MolFromSmarts("C(O)C(N)CO")
    # Pattern 2: A carbonyl variant (found in dehydro/unsaturated forms): C(=O)C(N)CO
    core_pattern2 = Chem.MolFromSmarts("C(=O)C(N)CO")
    
    # Now, to reduce the chance that the match is in a ring system (leading to false positives),
    # we require that at least one match (for either pattern) is entirely acyclic.
    def acyclic_match(pattern):
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # If none of the atoms in the match are in a ring, we consider this a valid hit.
            if all(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
                return True
        return False

    if not (acyclic_match(core_pattern1) or acyclic_match(core_pattern2)):
        return False, "No acyclic sphingoid core (amino-diol or carbonyl variant) found."
    
    # Check molecular weight: Many sphingoid molecules lie roughly between 200 to 800 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical sphingoid compound."
    
    # If all basic conditions are met, classify as sphingoid.
    return True, "Molecule exhibits an acyclic sphingoid core with a long aliphatic chain and appropriate functionality."

# Example test cases (for manual testing)
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCC\\C=C\\C(=O)[C@@H](N)CO",  # 3-dehydrosphingosine (should return True)
        "OC[C@@]([C@@](CCCCCCCCCCCCCC)(O)H)(N)H",  # Heptadecasphinganine (True)
        "CCO",  # too short (False)
        "CCCCCCCCCCC(=O)[C@@H](N)CO",  # 3-dehydrotetradecasphinganine (False under previous def., now may pass if other features match)
    ]
    for smi in test_smiles:
        classified, explanation = is_sphingoid(smi)
        print(f"SMILES: {smi}\n Classified: {classified}\n Reason: {explanation}\n")