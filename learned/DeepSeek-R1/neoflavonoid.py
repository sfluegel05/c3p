"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid (CHEBI:XXXXX)
"""
from rdkit import Chem
from rdkit.Chem import Mol, GetMolFrags, RWMol
from typing import Tuple

def is_neoflavonoid(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    benzopyran_found = False
    benzene_ring = None
    pyran_ring = None
    
    # Search for benzopyran core (benzene fused to a six-membered ring with one oxygen)
    for benzene in rings:
        # Check if benzene is a 6-membered aromatic carbon ring
        if len(benzene) != 6:
            continue
        is_aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in benzene)
        all_carbon = all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in benzene)
        if not (is_aromatic and all_carbon):
            continue
        
        # Look for fused six-membered ring with exactly one oxygen
        for other in rings:
            if other == benzene:
                continue
            if len(other) != 6:
                continue
            o_count = sum(1 for a in other if mol.GetAtomWithIdx(a).GetAtomicNum() == 8)
            if o_count != 1:
                continue
            
            # Check if the two rings share at least two adjacent atoms (fusion)
            shared = list(set(benzene) & set(other))
            if len(shared) < 2:
                continue
            
            # Check adjacency in benzene ring (any two consecutive atoms in benzene sharing with other ring)
            # Generate pairs of consecutive atoms in benzene
            benzene_pairs = [(benzene[i], benzene[(i+1)%6]) for i in range(6)]
            shared_pairs = [(shared[i], shared[j]) for i in range(len(shared)) for j in range(i+1, len(shared))]
            
            # Check if any shared pair is consecutive in benzene
            fused = any(pair in benzene_pairs for pair in shared_pairs)
            if fused:
                benzopyran_found = True
                benzene_ring = benzene
                pyran_ring = other
                break
        if benzopyran_found:
            break
    
    if not benzopyran_found:
        return False, "No 1-benzopyran core found"
    
    # Check for aryl substituent on pyran ring (not part of benzene)
    for atom_idx in pyran_ring:
        if atom_idx in benzene_ring:
            continue  # Skip fused atoms
        
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check all bonds from this atom
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetIdx() in pyran_ring:
                continue  # Part of the pyran ring
            
            # Check if the substituent contains an aromatic ring
            modified_mol = RWMol(mol)
            modified_mol.RemoveBond(atom.GetIdx(), neighbor.GetIdx())
            frags = GetMolFrags(modified_mol, asMols=True, sanitizeFrags=False)
            
            for frag in frags:
                if neighbor.GetIdx() not in [a.GetIdx() for a in frag.GetAtoms()]:
                    continue
                try:
                    Chem.SanitizeMol(frag)
                except:
                    continue
                # Check for aromatic atoms in the fragment
                if any(a.GetIsAromatic() for a in frag.GetAtoms()):
                    return True, "1-benzopyran with aryl substituent at position 4"
    
    return False, "No aryl substituent at position 4"