"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated fatty acyl-CoA
Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    That is, the molecule must (a) contain a thioester linkage typical of a fatty acyl-CoA,
    (b) contain a CoA fragment, and (c) its fatty acyl chain (the chain connected to the carbonyl
    carbon of the thioester) contains exactly one carbon-carbon double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check for a thioester group. Fatty acyl-CoAs have a thioester bond in the form of [C](=O)[S]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (lack of fatty acyl linkage)"
    
    # We assume the first thioester match corresponds to the linkage between the fatty acyl chain and CoA.
    # In a thioester match, atom0 is the carbonyl carbon.
    acyl_carbon_idx = thioester_matches[0][0]
    acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    
    # Check that the molecule has a CoA-related fragment.
    # Many CoA compounds share the fragment "SCCNC(=O)CCNC(=O)". This SMARTS is a heuristic.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Now, extract the fatty acyl chain atoms.
    # We start from the carbonyl carbon (which is part of the acyl chain) and perform a DFS 
    # but only follow bonds to carbon atoms (atomic number 6). This avoids traversing into the
    # thioester sulfur or other CoA portions.
    fatty_chain_atoms = set()
    stack = [acyl_carbon]
    while stack:
        atom = stack.pop()
        if atom.GetIdx() in fatty_chain_atoms:
            continue
        # only consider carbon atoms (atomic num 6)
        if atom.GetAtomicNum() != 6:
            continue
        fatty_chain_atoms.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            # Only traverse if neighbor is a carbon.
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in fatty_chain_atoms:
                stack.append(neighbor)
    
    if not fatty_chain_atoms:
        return False, "No fatty acyl chain detected"

    # Count carbon-carbon double bonds within the fatty chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check if both ends are in our fatty chain and the bond type is DOUBLE.
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        if a1 in fatty_chain_atoms and a2 in fatty_chain_atoms:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} carbon-carbon double bond(s), but exactly one is required"
    
    return True, "Contains a fatty acyl chain with exactly one C=C double bond and a CoA moiety"

# The function can be tested with one or more of the sample SMILES strings provided.
if __name__ == "__main__":
    # Example test with one of the provided SMILES (oleoyl-CoA)
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
    result, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)