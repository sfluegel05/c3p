"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: A polyunsaturated fatty acid that contains exactly three carbon–carbon double bonds,
with a terminal carboxylic acid group and an acyclic, long aliphatic chain.
This program uses several heuristics:
  1. The molecule must be parsed successfully (a valid structure).
  2. It must contain exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H1]").
  3. The carboxylic acid group must be terminal (its carbon is attached to only one carbon).
  4. There must be exactly three C=C (carbon–carbon) double bonds (ignoring C=O).
  5. There must be a sufficient number of carbon atoms (e.g. at least 10) and low cyclic content.
  6. The molecular weight should be above a minimal threshold.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    It verifies that the molecule:
      - Is valid.
      - Contains exactly one terminal carboxylic acid group.
      - Has exactly three carbon–carbon double bonds (excluding C=O bonds).
      - Contains enough carbons in an acyclic chain (low ring-fraction).
      - Has a molecular weight above a defined cutoff.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a trienoic fatty acid, False otherwise.
        str: A reason describing the classification result.
    """
    # Parse SMILES into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for exactly one carboxylic acid group.
    # SMARTS for carboxylic acid group: a carbonyl carbon double-bonded to an O and single bonded to an -OH.
    ca_smarts = "[CX3](=O)[OX2H1]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if len(ca_matches) == 0:
        return False, "Missing carboxylic acid functional group (COOH)"
    elif len(ca_matches) > 1:
        return False, f"Found {len(ca_matches)} carboxylic acid groups; requires exactly 1"
    
    # 2. Ensure the carboxylic acid group is terminal.
    # In a free fatty acid, the acid carbon should be bonded to only one carbon.
    ca_indices = ca_matches[0]
    ca_carbon_idx = ca_indices[0]  # Index of the carboxyl carbon in the match.
    ca_carbon = mol.GetAtomWithIdx(ca_carbon_idx)
    carbon_neighbors = [nbr for nbr in ca_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (expected to have exactly one carbon neighbor)"
    
    # 3. Count the number of carbon–carbon double bonds.
    # We count only those double bonds where both atoms are carbons,
    # and we exclude any double bond that may be directly part of the COOH group.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Exclude a bond involving the carboxyl carbon (since C=O here belongs to COOH)
                if ca_carbon_idx in (a1.GetIdx(), a2.GetIdx()):
                    continue
                double_bond_count += 1
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon–carbon double bonds; requires exactly 3 for a trienoic fatty acid"
    
    # 4. Check that the fatty acid chain is sufficiently long, based on number of carbons.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    total_carbons = len(carbons)
    if total_carbons < 10:
        return False, f"Too few carbon atoms ({total_carbons}); must be a long-chain fatty acid"
    
    # 5. Ensure that the molecule is largely acyclic.
    # We compute the fraction of carbon atoms that are part of any ring.
    ring_info = mol.GetRingInfo()
    ring_carbons = set()
    for ring in ring_info.AtomRings():
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6:
                ring_carbons.add(idx)
    fraction_ring = len(ring_carbons) / total_carbons if total_carbons else 0
    if fraction_ring > 0.2:
        return False, f"High cyclic content ({fraction_ring*100:.1f}% of carbons in rings); not a typical fatty acid chain"
    
    # 6. Check the molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180:
        return False, f"Molecular weight {mol_wt:.1f} Da too low for a fatty acid"
    
    return True, "Contains one terminal carboxylic acid group, exactly 3 carbon–carbon double bonds, and an acyclic long-chain structure typical of trienoic fatty acids"

# Example usage:
if __name__ == '__main__':
    # Test with one known trienoic fatty acid: 10,12,14-octadecatrienoic acid.
    test_smiles = "OC(=O)CCCCCCCC/C=C/C=C/C=C/CCC"
    result, reason = is_trienoic_fatty_acid(test_smiles)
    print(result, reason)