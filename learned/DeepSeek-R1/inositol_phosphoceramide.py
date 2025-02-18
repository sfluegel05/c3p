"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: CHEBI: inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    Must have inositol linked via phosphodiester to a ceramide with amide and long chains.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Generalized inositol pattern: cyclohexane with >=5 hydroxyls (including possible substitutions)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([OX2])[C@H]([OX2])[C@H]([OX2])[C@H]([OX2])[C@@H]([OX2])[C@@H]1[OX2]")
    if not mol.HasSubstructMatch(inositol_pattern):
        inositol_alt = Chem.MolFromSmarts("C1C(C(C(C(C1O)O)O)O)O")  # Less stereospecific
        if not mol.HasSubstructMatch(inositol_alt):
            return False, "No inositol-like structure"

    # Find phosphate connected to inositol
    phosphate = Chem.MolFromSmarts("[O]P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if not phosphate_matches:
        return False, "No phosphate group"

    # Check phosphodiester bridge connects inositol to ceramide
    inositol_atoms = set(mol.GetSubstructMatch(inositol_pattern) if mol.HasSubstructMatch(inositol_pattern) else mol.GetSubstructMatch(inositol_alt))
    bridge_found = False
    for p_match in phosphate_matches:
        p_idx = p_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        # Check phosphate has two ester linkages (O connected to inositol and another O to ceramide)
        ester_links = 0
        ceramide_link = False
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:  # Oxygen in ester
                # Check if connected to inositol
                for bond in neighbor.GetBonds():
                    other = bond.GetOtherAtom(neighbor)
                    if other.GetIdx() in inositol_atoms:
                        ester_links +=1
                    else:
                        # Follow this oxygen's connection to check for ceramide (amide with long chain)
                        stack = [(neighbor, 0)]
                        visited = set()
                        while stack:
                            current, depth = stack.pop()
                            if current.GetIdx() in visited or depth > 6:
                                continue
                            visited.add(current.GetIdx())
                            # Look for amide group (CONH)
                            if current.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in current.GetBonds()):
                                for nbr in current.GetNeighbors():
                                    if nbr.GetAtomicNum() == 7 and any(bond.GetBondType() == Chem.BondType.SINGLE for bond in nbr.GetBonds()):
                                        # Check for long chain (>=14 carbons)
                                        chain_length = 0
                                        chain_stack = [(current, 0)]
                                        chain_visited = set()
                                        while chain_stack:
                                            catom, clevel = chain_stack.pop()
                                            if catom.GetIdx() in chain_visited:
                                                continue
                                            chain_visited.add(catom.GetIdx())
                                            if catom.GetAtomicNum() == 6:
                                                chain_length +=1
                                                if chain_length >=14:
                                                    ceramide_link = True
                                                    break
                                                for b in catom.GetBonds():
                                                    next_a = b.GetOtherAtom(catom)
                                                    if next_a.GetAtomicNum() == 6 and clevel <30:
                                                        chain_stack.append((next_a, clevel+1))
                                            if ceramide_link:
                                                break
                            if not ceramide_link:
                                for bond in current.GetBonds():
                                    next_a = bond.GetOtherAtom(current)
                                    if next_a.GetAtomicNum() in [6,8]:
                                        stack.append((next_a, depth+1))
        if ester_links >=1 and ceramide_link:
            bridge_found = True
            break
    if not bridge_found:
        return False, "No phosphodiester bridge to ceramide"

    # Verify ceramide has sphingoid base with hydroxyl
    amide_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[NX3H]"))
    if not amide_matches:
        return False, "No amide group"
    hydroxyl_found = False
    for amide in amide_matches:
        nh_atom = mol.GetAtomWithIdx(amide[2])
        # Check adjacent to NH for hydroxyl (sphingoid base)
        for neighbor in nh_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >0:
                hydroxyl_found = True
                break
        if not hydroxyl_found:
            # Check within 3 bonds for hydroxyl
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, nh_atom.GetIdx())
            atoms = set()
            for bidx in env:
                atoms.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())
            for aidx in atoms:
                a = mol.GetAtomWithIdx(aidx)
                if a.GetAtomicNum() ==8 and a.GetTotalNumHs()>0:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "No hydroxyl in sphingoid base"

    return True, "Inositol-phosphoceramide with phosphodiester bridge and ceramide features"