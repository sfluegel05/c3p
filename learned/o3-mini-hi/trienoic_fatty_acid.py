"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three carbon‐carbon double bonds.
A genuine trienoic fatty acid must have one (and only one) carboxylic acid group,
exactly three C=C bonds between carbon atoms (excluding C=O bonds in the acid),
a sufficient number of carbon atoms (indicating a long, aliphatic chain),
and a low fraction of carbons in rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    The method applies several heuristics:
      1. The molecule must be a valid structure.
      2. It has exactly one carboxylic acid group (matched by SMARTS "[CX3](=O)[OX2H1]").
      3. It has exactly three carbon–carbon double bonds (excluding C=O bonds).
      4. It has a minimum number of carbon atoms (heuristically at least 8).
      5. It has a low fraction of carbon atoms involved in rings (indicating an acyclic, aliphatic chain).
         (Here we use a threshold of 20% of all carbons.)
      6. Additionally, we check that the molecular weight is above a minimal threshold.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a trienoic fatty acid, False otherwise.
        str: A reason describing the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for exactly one carboxylic acid group.
    # SMARTS: a carbon with 3 connections, double bonded to oxygen and connected to -OH.
    ca_smarts = "[CX3](=O)[OX2H1]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if len(ca_matches) == 0:
        return False, "Missing carboxylic acid functional group (COOH)"
    elif len(ca_matches) > 1:
        return False, f"Found {len(ca_matches)} carboxylic acid groups; requires exactly 1 for a fatty acid"
    
    # 2. Count the number of carbon-carbon double bonds.
    # Only count bonds that are truly C=C (i.e. both atoms are carbon).
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Skip double bonds where one atom is not carbon (e.g. C=O in the acid)
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon-carbon double bonds; requires exactly 3 for trienoic fatty acid"
    
    # 3. Ensure a minimal number of carbon atoms are present.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    total_carbons = len(carbon_atoms)
    if total_carbons < 8:
        return False, f"Too few carbon atoms ({total_carbons}); must be a long-chain fatty acid"
    
    # 4. Check the fraction of carbon atoms that are in rings.
    ring_info = mol.GetRingInfo()
    # Gather the set of carbon atoms that are in any ring.
    ring_carbon_ids = set()
    for ring in ring_info.AtomRings():
        for atom_id in ring:
            atom = mol.GetAtomWithIdx(atom_id)
            if atom.GetAtomicNum() == 6:
                ring_carbon_ids.add(atom_id)
    fraction_in_rings = len(ring_carbon_ids) / total_carbons if total_carbons > 0 else 0
    # Fatty acids are expected to be mostly acyclic.
    if fraction_in_rings > 0.2:
        return False, f"High degree of cyclic (ring) carbons ({fraction_in_rings*100:.1f}%); not a typical fatty acid chain"
    
    # 5. Check molecular weight (heuristic: fatty acids are not extremely light).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da too low for a fatty acid"
    
    return True, "Contains one carboxylic acid group, exactly 3 carbon-carbon double bonds, and an acyclic long-chain structure typical of trienoic fatty acids"

# Example usage:
if __name__ == '__main__':
    test_smiles = "OC(=O)CCCCCCCC/C=C/C=C/C=C/CCC"  # 10,12,14-octadecatrienoic acid
    result, reason = is_trienoic_fatty_acid(test_smiles)
    print(result, reason)