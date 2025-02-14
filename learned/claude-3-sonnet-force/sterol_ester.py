"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:49487 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group
    of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid backbone patterns
    steroid_patterns = [Chem.MolFromSmarts("[C@@]12[C@H](CC[C@@]1(C)CCC3=C2CC[C@H]4[C@@]3(CCC(C4)C)C)"),
                        Chem.MolFromSmarts("[C@@]12[C@H](CC[C@@]1(C)CCC3=CC(=O)[C@H]4[C@@]3(CCC(C4)C)C)"),
                        Chem.MolFromSmarts("[C@@]12[C@H](CC[C@]1([H])CCC3=C2CC[C@@]4([H])[C@]3(CCC(C4)C)C)")]
    
    steroid_match = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(pattern):
            steroid_match = True
            break
            
    if not steroid_match:
        return False, "No steroid backbone found"
    
    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check if ester group is connected to steroid backbone
    steroid_atoms = set()
    for match in ester_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        if any(mol.GetAtomWithIdx(nb).GetAtomicNum() == 6 and nb not in match for nb in o_atom.GetNeighbors()):
            steroid_atoms.update(o_atom.GetNeighbors())
    
    # Look for fatty acid chain (long carbon chain attached to ester)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = []
    for match in ester_matches:
        c_atom = mol.GetAtomWithIdx(match[1])
        for nb in c_atom.GetNeighbors():
            if nb not in match and mol.GetAtomWithIdx(nb).GetAtomicNum() == 6:
                chain_match = mol.GetSubstructMatches(chain_pattern, atomIds=[nb])
                if chain_match:
                    chain_matches.extend(chain_match)
                    
    if not chain_matches:
        return False, "No fatty acid chain found"
    
    # Check if fatty acid chain is connected to steroid backbone
    chain_atoms = set([m[0] for m in chain_matches])
    if not chain_atoms.intersection(steroid_atoms):
        return False, "Fatty acid chain not connected to steroid backbone"
    
    # Count rotatable bonds and molecular weight
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    if n_rotatable < 8 or mol_wt < 400:
        return False, "Molecule too small to be a sterol ester"
    
    return True, "Contains steroid backbone with a fatty acid chain attached via an ester group"