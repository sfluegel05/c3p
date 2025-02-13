"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycoside group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate all possible tautomers and resonance structures
    tautomers = list(AllChem.ResonanceEnumerator(mol))
    tautomers.append(mol)

    # Enumerate all possible stereoisomers of D-glucose
    d_glucose_smarts = "[C@H]1([C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O"
    d_glucose_mols = [Chem.MolFromSmarts(d_glucose_smarts)]
    for perm in AllChem.EnumerateStereoisomers(Chem.MolFromSmarts(d_glucose_smarts)):
        d_glucose_mols.append(Chem.MolFromMolBlock(Chem.MolToMolBlock(perm)))

    for tautomer in tautomers:
        for d_glucose in d_glucose_mols:
            # Check if tautomer contains the D-glucose substructure
            if tautomer.HasSubstructMatch(d_glucose):
                # Find glycosidic bonds involving the D-glucose substructure
                for match in tautomer.GetSubstructMatches(d_glucose):
                    d_glucose_atoms = list(match)
                    for atom_idx in d_glucose_atoms:
                        atom = tautomer.GetAtomWithIdx(atom_idx)
                        for bond in atom.GetBonds():
                            if bond.GetBondType() == BondType.SINGLE and bond.GetOtherAtomIdx(atom_idx) not in d_glucose_atoms:
                                neighbor_atom = tautomer.GetAtomWithIdx(bond.GetOtherAtomIdx(atom_idx))
                                if neighbor_atom.GetSymbol() == "O" and neighbor_atom.GetTotalNumHs() == 0:
                                    for bond2 in neighbor_atom.GetBonds():
                                        if bond2.GetBondType() == BondType.SINGLE and bond2.GetOtherAtomIdx(neighbor_atom.GetIdx()) not in d_glucose_atoms:
                                            # Found a glycosidic bond involving D-glucose
                                            return True, "Contains a D-glucose substructure with a glycosidic bond"

    return False, "No D-glucose substructure with a glycosidic bond found"