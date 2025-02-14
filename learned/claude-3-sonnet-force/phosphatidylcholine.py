"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: CHEBI:17556 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine with two acyl substituents at positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
        glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
        if not mol.HasSubstructMatch(glycerol_pattern):
            return False, "No glycerol backbone found"

        # Look for phosphocholine group (-OP(=O)([O-])OCC[N+](C)(C)C)
        phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
        if not mol.HasSubstructMatch(phosphocholine_pattern):
            return False, "No phosphocholine group found"

        # Identify ester groups
        ester_matches = Chem.AddHs(mol).GetSubstructMatches(Chem.MolFromSmarts('COC=O'))
        if len(ester_matches) != 2:
            return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

        # Check for long carbon chains (fatty acid chains) attached to esters
        fatty_acid_chains = []
        for ester_idx in ester_matches:
            ester_atom = mol.GetAtomWithIdx(ester_idx)
            for neighbor in ester_atom.GetNeighbors():
                chain = []
                current_atom = neighbor
                while current_atom.GetAtomicNum() == 6:
                    chain.append(current_atom.GetIdx())
                    neighbors = current_atom.GetNeighbors()
                    if len(neighbors) == 1:
                        break
                    current_atom = [n for n in neighbors if n.GetIdx() != current_atom.GetPrevAtomIdx()][0]
                if len(chain) >= 4:
                    fatty_acid_chains.append(chain)

        if len(fatty_acid_chains) != 2:
            return False, "Missing fatty acid chains"

        # Count rotatable bonds to verify long chains
        n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if n_rotatable < 10:
            return False, "Chains too short to be fatty acids"

        # Calculate expected molecular weight based on detected fatty acid chains
        expected_mol_wt = rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles('C(CCCCCCCCCCCCCCCCCCCCC)OP(=O)([O-])OCC[N+](C)(C)C'))
        for chain in fatty_acid_chains:
            for idx in chain:
                atom = mol.GetAtomWithIdx(idx)
                expected_mol_wt += atom.GetMass() - 12.011  # Subtract carbon mass

        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if abs(mol_wt - expected_mol_wt) > 10:
            return False, "Molecular weight does not match expected value for phosphatidylcholine"

        return True, "Contains glycerol backbone with 2 fatty acid chains and phosphocholine group"

    except Exception as e:
        return False, f"Error: {str(e)}"