"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
The function identifies the amino acid backbone by searching for an α‐carbon connected to both 
an amino group and a carboxyl group. It then extracts the side chain and looks for polar groups 
(i.e. those containing O or N, or an –SH group with at least one hydrogen).
"""

from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines whether a given amino acid (as a SMILES string) is a polar amino acid, meaning
    that its side chain has at least one functional group capable of hydrogen bonding.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a polar amino acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Our goal is to identify the amino acid backbone.
    # We look for an α-carbon (atom C) that is connected to:
    #   1) An amino group (N) and
    #   2) A carboxyl group (a carbon that is double-bonded to an oxygen, i.e. C=O).
    alpha_idx = None          # index of the α-carbon
    sidechain_start_idx = None  # the neighbor (atom index) representing the start of the side chain

    for atom in mol.GetAtoms():
        # We consider only carbon atoms as possible α-carbons.
        if atom.GetAtomicNum() != 6:
            continue
        # In an amino acid with implicit hydrogens, the heavy-atom degree of the α-carbon is usually 3.
        if len(atom.GetNeighbors()) != 3:
            continue
        
        neighbors = atom.GetNeighbors()
        has_amino = False
        carboxyl_idx = None
        
        # Identify an amino neighbor: a nitrogen atom.
        for nb in neighbors:
            if nb.GetAtomicNum() == 7:
                has_amino = True
                break
        if not has_amino:
            continue
        
        # Identify an attached carboxyl carbon.
        for nb in neighbors:
            # Looking for a carbon, not the α-carbon itself, that is part of a carboxyl group.
            if nb.GetAtomicNum() == 6:
                # Check if this neighbor has an oxygen with a double bond.
                found_carboxyl = False
                for bond in nb.GetBonds():
                    # Get the other atom in the bond.
                    other = bond.GetOtherAtom(nb)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        found_carboxyl = True
                        break
                if found_carboxyl:
                    carboxyl_idx = nb.GetIdx()
                    break
        if carboxyl_idx is None:
            continue
        
        # Having found an α-carbon candidate having at least one amino neighbor and one carboxyl neighbor,
        # the remaining neighbor is considered as the side chain starting point.
        for nb in neighbors:
            if nb.GetAtomicNum() == 7:
                continue
            if nb.GetIdx() == carboxyl_idx:
                continue
            # This neighbor is the start of the side chain.
            sidechain_start_idx = nb.GetIdx()
            break
            
        if sidechain_start_idx is not None:
            # We have identified the backbone; record α-carbon index.
            alpha_idx = atom.GetIdx()
            break

    if alpha_idx is None or sidechain_start_idx is None:
        return False, "Amino acid backbone (α-carbon with amino and carboxyl groups) not detected"
    
    # Extract the side chain substructure by performing a DFS starting from the side chain atom,
    # ensuring we do not traverse back into the backbone (α-carbon).
    sidechain_atoms = set()
    stack = [sidechain_start_idx]
    while stack:
        current_idx = stack.pop()
        if current_idx in sidechain_atoms:
            continue
        sidechain_atoms.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nb in current_atom.GetNeighbors():
            # Do not traverse back into the backbone (skip the α-carbon atom)
            if nb.GetIdx() == alpha_idx:
                continue
            if nb.GetIdx() not in sidechain_atoms:
                stack.append(nb.GetIdx())
    
    # Now analyze the side chain fragment for the presence of polar atoms.
    # We treat any oxygen (atomic num 8) or nitrogen (atomic num 7) as polar. 
    # In addition, for sulfur (atomic num 16) we require that at least one hydrogen is attached (i.e. –SH).
    polar_found = False
    polar_features = []
    for idx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        anum = atom.GetAtomicNum()
        if anum == 7 or anum == 8:
            # Nitrogen or oxygen in the side chain: consider polar.
            polar_found = True
            polar_features.append(f"{atom.GetSymbol()}")
        elif anum == 16:
            # For sulfur, check that at least one hydrogen is attached.
            if atom.GetTotalNumHs() > 0:
                polar_found = True
                polar_features.append("SH")
    
    if polar_found:
        return True, f"Side chain contains polar feature(s): {', '.join(polar_features)}"
    else:
        return False, "Side chain does not contain a polar functional group capable of hydrogen bonding"

# Example usage (this portion can be omitted or used for testing):
if __name__ == "__main__":
    test_smiles = [
        "NC(CCC(N)=O)C(O)=O",   # glutamine, polar
        "NC(CO)C(O)=O",        # serine, polar
        "NCC(=O)O",            # glycine (side chain = H, non-polar)
        "NC(CS)C(O)=O",        # cysteine, polar due to –SH
    ]
    
    for s in test_smiles:
        result, reason = is_polar_amino_acid(s)
        print(f"SMILES: {s}\nIs polar amino acid? {result} ({reason})\n")