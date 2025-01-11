"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid 
with a hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a tetradecanoate ester group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for non-cyclic ester group pattern
    # [C:1] must be sp3 carbon (no double bonds)
    # Must not be in ring
    ester_pattern = Chem.MolFromSmarts("[CH3X4:1][CH2X4:2][CH2X4:3][CH2X4:4][CH2X4:5][CH2X4:6][CH2X4:7][CH2X4:8][CH2X4:9][CH2X4:10][CH2X4:11][CH2X4:12][CH2X4:13][CX3:14](=[OX1:15])[OX2:16][*:17]")
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No tetradecanoate ester group found"
        
    matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in matches:
        # Get the matched atoms
        chain_atoms = match[0:13]  # First 13 carbons
        carbonyl_carbon = match[13]
        ester_oxygen = match[15]
        
        # Verify none of the chain atoms are in a ring
        ring_atoms = mol.GetRingInfo().AtomRings()
        if any(atom_idx in {atom for ring in ring_atoms for atom in ring} for atom_idx in chain_atoms):
            continue
            
        # Verify chain carbons don't have additional carbon substituents
        has_branches = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
            if carbon_neighbors > 2:  # More than 2 carbon neighbors means branching
                has_branches = True
                break
                
        if has_branches:
            continue
            
        # Verify no modifications on the chain (OH groups, double bonds, etc)
        has_modifications = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetTotalNumHs() + len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) < 4:
                has_modifications = True
                break
                
        if has_modifications:
            continue
            
        # If we get here, we've found a valid tetradecanoate ester group
        return True, "Contains tetradecanoate (myristoyl) ester group"
        
    return False, "No valid tetradecanoate ester group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:87691',
                          'name': 'tetradecanoate ester',
                          'definition': 'A fatty acid ester obtained by '
                                        'condensation of the carboxy group of '
                                        'tetradecanoic acid (also known as '
                                        'myristic acid) with a hydroxy group '
                                        'of an alcohol or phenol.',
                          'parents': ['CHEBI:35748'],
                          'xrefs': ['PMID:26212120'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.48044692737430167 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\\C=C/CCCCCCCC)COP(O)(O)=O '
               'NAME: 1-oleoyl-2-myristoyl-sn-glycero-3-phosphate REASON: '
               'CORRECT Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@@H](COC(=O)CCCCCCC\\C=C/CCCCCCCC)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1,2-dioleoyl-3-myristoyl-sn-glycerol REASON: CORRECT '
               'Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C '
               'NAME: O-tetradecanoyl-L-carnitine REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-myristoyl-sn-glycerol REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@@H](O)CO NAME: '
               "1-tetradecanoyl-sn-glycero-3-phospho-(1'-sn-glycerol) REASON: "
               'CORRECT Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC '
               'NAME: tetramyristoyl cardiolipin REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@@H](CO)COC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-oleoyl-2-myristoyl-sn-glycerol REASON: CORRECT '
               'Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-myristoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine '
               'REASON: CORRECT Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C NAME: '
               'O-tetradecanoylcarnitine REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[S+](C)C)OC(=O)CCCCCCCCCCCCC '
               'NAME: 1,2-di-O-myristoyl-sn-glycero-3-phosphosulfocholine '
               'REASON: CORRECT Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: C(CCCCCCCC)CCCCC(OCCCCCCCCCCCCCC)=O NAME: '
               'tetradecyl tetradecanoate REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\\C=C/CCCCCCCC)COP([O-])(=O)OCC[NH3+] '
               'NAME: 1-oleoyl-2-myristoyl-sn-glycero-3-phosphoethanolamine '
               'zwitterion REASON: CORRECT Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](CO)OC(=O)CCCCCCC\\C=C/CCCCCC NAME: '
               '1-myristoyl-2-palmitoleoyl-sn-glycerol REASON: CORRECT '
               'Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](CCCCCCCCCCC)CC(=O)O[C@@H]1[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@H](OC[C@H]2O[C@H](OP(O)(O)=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO[C@@]2(C[C@@H](O[C@@]3(C[C@@H](O)[C@@H](O)[C@H](O3)[C@H](O)CO)C(O)=O)[C@@H](O)[C@H](O2)[C@H](O)CO)C(O)=O)[C@H]1OP(O)(O)=O '
               'NAME: (Kdo)2-lipid A (E. coli) REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC '
               'NAME: '
               '1-tetradecanoyl-2-[(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl]-sn-glycero-3-phosphocholine '
               'REASON: CORRECT Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCCCCCCCC NAME: '
               '1,2 ditetradecanoylglycerol REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC '
               'NAME: 1,2-di-O-myristoyl-sn-glycero-3-phosphocholine REASON: '
               'CORRECT Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](CCCCCCCCCCC)CC(=O)O[C@@H]1[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@H](OC[C@H]2O[C@H](OP(O)(=O)OP(O)(=O)OCCN)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO)[C@H]1OP(O)(=O)OP(O)(=O)OCCN '
               "NAME: lipid A 1,4'-bis(2-aminoethyl diphosphate) REASON: "
               'CORRECT Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCC '
               'NAME: 1-myristoyl-2-palmitoleoyl-sn-glycero-3-phosphocholine '
               'REASON: CORRECT Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H]1CC[C@@](C)(OC(C)=O)[C@@H]2O[C@H](C[C@]1(C)O)[C@@H]1[C@H]2[C@@H](C(C)C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@@]11CO1 '
               'NAME: klymollin F REASON: CORRECT Contains tetradecanoate '
               '(myristoyl) ester group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC[C@@H](CO)OC(=O)CCCCCCCCCCCCC '
               'NAME: 2,3-dimyristoyl-sn-glycerol REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC NAME: '
               '1-myristoyl-2-oleoyl-sn-glycerol REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[NH3+])OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: 1-myristoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine '
               'zwitterion REASON: CORRECT Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC[C@@H](O)COP([O-])(=O)OCC[NH3+] '
               'NAME: 1-myristoyl-sn-glycero-3-phosphoethanolamine zwitterion '
               'REASON: CORRECT Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC NAME: '
               '1-myristoyl-2-oleoylglycerol REASON: CORRECT Contains '
               'tetradecanoate (myristoyl) ester group\n'
               'False positives: SMILES: '
               'CCCCCCCCCCC[C@@H](O)CC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@H](OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@@H]1NC(C)=O '
               'NAME: '
               'UDP-3-O-[(3R)-3-hydroxytetradecanoyl]-N-acetylglucosamine(2-) '
               'REASON: WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: O=C1O[C@H](C[C@H]([C@@H]1O)O)CCCCCCCCC NAME: '
               '(2S,3R,5S)-(-)-2,3-dihydroxytetradecan-5-olide REASON: WRONGLY '
               'CLASSIFIED Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C NAME: '
               '(5Z,8Z)-tetradecadienoylcarnitine REASON: WRONGLY CLASSIFIED '
               'Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@@H](OC(=O)/C(=C/C(CC)C)/C)[C@@H](O)[C@@H]2[C@@H]([C@H]([C@@H](C1)O)O)O2)CCCCC '
               'NAME: 5,6-Epoxy-phomol REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CC(O)CC(=O)OC(CC([O-])=O)C[N+](C)(C)C NAME: '
               'O-[(5Z,8Z)-3-hydroxytetradecadienoyl]carnitine REASON: WRONGLY '
               'CLASSIFIED Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: O=C1O[C@@H]([C@@H](O)/C=C/C=C\\CC[C@H](O)CC)CC1 '
               'NAME: Butanolide A REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)O[C@H](CCCCCCCCCCC)CC(=O)O[C@@H]1[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@H](OC[C@H]2O[C@H](OP(O)(O)=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO[C@@]2(C[C@@H](O[C@@]3(C[C@@H](O)[C@@H](O)[C@H](O3)[C@H](O)CO)C(O)=O)[C@@H](O[C@H]3O[C@H]([C@H](CO)OP(O)(=O)OCCCN)[C@@H](O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)[C@H](O[C@H]4O[C@H]([C@H](O)CO)[C@@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@H](O[C@@H]6O[C@H](CO)[C@H](O)[C@H](O[C@@H]7O[C@H](CO)[C@H](O)[C@H](O[C@@]8(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O8)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]7O)[C@H]6NC(C)=O)[C@H](O[C@@]6(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)[C@@H]4O)[C@@H]3O)[C@H](O2)[C@H](O)CO)C(O)=O)[C@H]1OP(O)(O)=O '
               'NAME: '
               'alpha-Neup5Ac-(2->3)-beta-D-Galp-(1->3)-beta-D-GalpNAc-(1->4)-[alpha-Neup5Ac-(2->3)]-beta-D-Galp-(1->3)-L-alpha-D-Hepp-(1->3)-[beta-D-Glcp-(1->4)]-L-alpha-D-Hepp6PEtn-(1->5)-[alpha-D-Kdo-(2->4)]-alpha-D-Kdo-(2->6)-lipid '
               'A REASON: WRONGLY CLASSIFIED Contains tetradecanoate '
               '(myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCC(O)CC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@@H](OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@@H]1NC(C)=O '
               'NAME: '
               'UDP-3-O-(3-hydroxytetradecanoyl)-N-acetyl-beta-glucosamine '
               'REASON: WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCC(O)CC(=O)O[C@H]1[C@H](O)[C@@H](CO)O[C@@H](OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@@H]1NC(C)=O '
               'NAME: '
               'UDP-3-O-(3-hydroxytetradecanoyl)-N-acetyl-beta-glucosamine(2-) '
               'REASON: WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: O=C1O[C@H]([C@H](O)[C@H](O)/C=C\\CCCCC(O)C)CC1 '
               'NAME: Pseudonocardide B REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H]1O[C@H](OC2CC[C@]3(C)[C@H]4CC[C@]5(C)[C@H](CC[C@H]5[C@@H]4CC=C3C2)[C@H](C)CCCC(C)C)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: cholesteryl 6-O-myristoyl-alpha-D-glucoside REASON: '
               'WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: O=C1OC(CC=C[C@H](CCC1)O)CCCCC NAME: Putaminoxin B '
               'REASON: WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](CCCCCCCCCCC)CC(=O)O[C@@H]1[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@H](OC[C@H]2O[C@H](OP([O-])(=O)OP([O-])([O-])=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO[C@@]2(C[C@@H](O[C@@]3(C[C@@H](O)[C@@H](O)[C@H](O3)[C@H](O)CO)C([O-])=O)[C@@H](O)[C@H](O2)[C@H](O)CO)C([O-])=O)[C@H]1OP([O-])([O-])=O '
               'NAME: (Kdo)2-lipid A 1-diphosphate(7-) (E. coli) REASON: '
               'WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC '
               'NAME: dimyristoyl phosphatidylcholine REASON: WRONGLY '
               'CLASSIFIED Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: C(CCCCCCCC)CCCCC(OC)=O NAME: methyl tetradecanoate '
               'REASON: WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) '
               'ester group\n'
               ' * SMILES: O=C1O[C@H]([C@H](O)/C=C/C(O)CCCCCC)CC1 NAME: '
               'Pseudonocardide E REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])([O-])=O)OC(=O)CCCCCCCCCCCCC '
               'NAME: 1,2-ditetradecanoyl-sn-glycerol-3-phosphate(2-) REASON: '
               'WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\\C=C/CCCCCCCC)COP([O-])([O-])=O '
               'NAME: 1-oleoyl-2-myristoyl-sn-glycero-3-phosphate(2-) REASON: '
               'WRONGLY CLASSIFIED Contains tetradecanoate (myristoyl) ester '
               'group\n'
               ' * SMILES: O=C1O[C@H]([C@H](O)[C@H](O)/C=C\\CCCCCC)[C@H](C1)O '
               'NAME: Pseudonocardide C REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: lysophosphatidylcholine 14:1(9Z)/0:0 REASON: WRONGLY '
               'CLASSIFIED Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@H](O)[C@@]([H])(CC(=O)O[C@@H](C)CCC\\C=C\\1)O2 '
               'NAME: aspergillide A REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: O=C1O[C@H](CCC[C@@H](OC)C=CC=CC(CC1)=O)C NAME: '
               'Pestalotioprolide H REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OCC(COP([O-])([O-])=O)OC(=O)CCCCCCCCCCCCC '
               'NAME: 1,2-dimyrsitoylphosphatidate(2-) REASON: WRONGLY '
               'CLASSIFIED Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)O[C@H](CCCCCCCCCCC)CC(=O)O[C@@H]1[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@H](OC[C@H]2O[C@H](OP([O-])([O-])=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO[C@@]2(C[C@@H](O[C@@]3(C[C@@H](O)[C@@H](O)[C@H](O3)[C@H](O)CO)C([O-])=O)[C@@H](O)[C@H](O2)[C@H](O)CO)C([O-])=O)[C@H]1OP([O-])([O-])=O '
               'NAME: (Kdo)2-lipid A(6-) (E. coli) REASON: WRONGLY CLASSIFIED '
               'Contains tetradecanoate (myristoyl) ester group\n'
               ' * SMILES: O=C1O[C@H](C#CC#CC#CC#CC(O)C)[C@H](C1)O NAME: '
               '3,4,13-trihydroxy-tetradeca-5,7,9,11-tetraynoic '
               'acid-gamma-lac-tone REASON: WRONGLY CLASSIFIED Contains '
               'tetradecanoate (myristoyl) ester group\n'
               'False negatives: SMILES: '
               '[C@H]1(OP(=O)(O)O)[C@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](O[C@@H]1CO)OC[C@@H]2[C@H]([C@@H]([C@H]([C@H](O2)OP(=O)(O)O)NC(=O)C[C@@H](CCCCCCCCCCC)OC(CCCCCCCCCCCCCCC)=O)OC(=O)C[C@@H](CCCCCCCCCCC)O)O '
               'NAME: palmitoyllipid A REASON: MISSED No tetradecanoate ester '
               'group found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCC)COP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: 1-tetradecanoyl-2-docosanoyl-sn-glycero-3-phosphocholine '
               'REASON: MISSED No tetradecanoate ester group found\n'
               ' * SMILES: '
               'C(C[N+](C)(C)C)OP(=O)([O-])OC[C@](COC(=O)CCCCCCCCCCCCC)([H])O '
               'NAME: 1-myristoyl-sn-glycero-3-phosphocholine REASON: MISSED '
               'No tetradecanoate ester group found\n'
               ' * SMILES: '
               '[C@H](OC(CCCCCCCCCCCCC)=O)(C[N+](C([2H])([2H])[2H])(C)C)CC(=O)[O-] '
               'NAME: tetradecanoyl-L-carnitine-d3 REASON: MISSED No '
               'tetradecanoate ester group found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC '
               'NAME: 1-palmitoyl-2-myristoyl-sn-glycero-3-phosphocholine '
               'REASON: MISSED No tetradecanoate ester group found\n'
               ' * SMILES: '
               '[C@H]1(OP(=O)(O)O)[C@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@H](OC[C@@H]2[C@H]([C@@H]([C@H]([C@H](O2)OP(=O)(O)O)NC(=O)C[C@@H](CCCCCCCCCCC)O)OC(=O)C[C@@H](CCCCCCCCCCC)O)O)O[C@@H]1CO[C@@]3(C[C@@H](O[C@]4(C(O)=O)C[C@@H](O)[C@H]([C@](O4)([C@H](O)CO)[H])O)[C@@H](O[C@@H]5[C@H]([C@H]([C@@H]([C@H](O5)[C@H](CO)O)O)O[C@@H]6[C@H]([C@H]([C@@H]([C@H](O6)[C@H](CO)O)O)O[C@@H]7[C@@H]([C@H]([C@@H]([C@H](O7)CO)O)O)O)O)O)[C@](O3)([C@@H](CO)O)[H])C(O)=O '
               'NAME: glucosyl-(heptosyl)2-(KDO)2-lipid A REASON: MISSED No '
               'tetradecanoate ester group found\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@@](C[C@@H](O[C@@]2(C[C@@H](O)[C@@H](O)[C@]([H])(O2)[C@@H](CO)OP(O)(=O)OCCN)C(O)=O)[C@H]1O)(OC[C@H]1O[C@@H](OC[C@H]2O[C@H](OP(O)(O)=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)[C@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H]1OP(O)(O)=O)C(O)=O)[C@H](O)CO '
               'NAME: phosphoethanolamine-Kdo2-lipid A REASON: MISSED No '
               'tetradecanoate ester group found\n'
               ' * SMILES: '
               '[C@H]1(OP(=O)(O)O[C@H]2OC[C@@H]([C@@H]([C@H]2O)O)N)[C@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](O[C@@H]1CO)OC[C@@H]3[C@H]([C@@H]([C@H]([C@H](O3)OP(=O)(O)O)NC(=O)C[C@@H](CCCCCCCCCCC)O)OC(=O)C[C@@H](CCCCCCCCCCC)O)O '
               'NAME: beta-L-Ara4N-lipid A (E. coli) REASON: MISSED No '
               'tetradecanoate ester group found\n'
               ' * SMILES: '
               'C=1[C@]2([C@]3([C@@]([C@H](OC(=O)CCCCCCCCCCCCC)[C@H]([C@@]2([C@]4([C@@](C(=O)C(=C4)C)(CC1CO)O)[H])O)C)(C3(C)C)OC(=O)C)[H])[H] '
               'NAME: phorbol 13-acetate 12-myristate REASON: MISSED No '
               'tetradecanoate ester group found\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)(O[C@@H]1C[C@@](OC[C@H]2O[C@@H](OC[C@H]3O[C@H](OP(O)(O)=O)[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]3O)[C@H](NC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](OC(=O)C[C@@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@@H]2OP(O)(O)=O)(O[C@]([H])([C@H](O)CO)[C@@H]1O[C@H]1O[C@H]([C@@H](O)CO)[C@@H](OP(O)(O)=O)[C@H](O[C@H]2O[C@H]([C@@H](O)COC3O[C@H]([C@@H](O)CO)[C@@H](O)[C@H](O)[C@@H]3O)[C@@H](OP(O)(O)=O)[C@H](OC3O[C@H](COC4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)[C@@H](O)[C@H](OC4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4OC4O[C@H](COC5O[C@H]([C@@H](O)CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O)[C@H]4O)[C@H]3O)[C@@H]2O)[C@@H]1O)C(O)=O)C(O)=O)[C@H](O)CO '
               'NAME: lipid A-core REASON: MISSED No tetradecanoate ester '
               'group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O1[C@H](CC2=C(C=3C=4C(C(O)=C5C3C[C@H](C)OC5)=C(OC)C=C(C4)OC)C6=CC(OC)=CC(=C6C(=C2C1)O)OC)C',
                                     'name': '(3S)-5-[(3S)-10-hydroxy-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-5-yl]-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-10-ol',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3CCC3',
                                     'name': 'N-[(5R,6S,9R)-8-[cyclobutyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': '[C@@]12(CCCC[C@@]1(C=C[C@@H]([C@@H]2CC[C@H](C[C@H](CC(=O)[O-])O)O)C)[H])[H]',
                                     'name': '4a,5-dihydro-ML-236C carboxylate',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'S(=O)(CC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1O)C',
                                     'name': "(S)-5'-Deoxy-5'-(methylsulfinyl)adenosine",
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'OC=1C(=C(C=2C=3C(NC2)=CC=CC3)C(=O)C(=O)C1C=4C=5C(NC4)=CC=CC5)CC=C(C)C',
                                     'name': 'Ochrindole D',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'O(C1C(O)C(OC1OC2=C(OC=3C(C2=O)=C(O)C=C(OC4OC(C(O)C(O)C4O)C)C3)C5=CC=C(O)C=C5)CO)C6OCC(O)(C6O)CO',
                                     'name': 'Kaempferol '
                                             '3-apiosyl-(1->2)-alpha-L-arabinofuranoside-7-rhamnoside',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': '[O-]C(=O)OON=O',
                                     'name': 'nitrosoperoxycarbonate(1-)',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'LysoPA(18:0/0:0)',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCC)CO/C=C\\CCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                     'name': 'PG(P-18:0/19:0)',
                                     'reason': 'No tetradecanoate ester group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Looking at the outcomes, there are several issues with the '
                 'previous code that need to be addressed:\n'
                 '\n'
                 '1. False Positives: The code is incorrectly identifying '
                 "molecules that contain 14-carbon chains but aren't "
                 'tetradecanoate esters, such as:\n'
                 '   - Lactones (cyclic esters)\n'
                 '   - Esters of hydroxy-tetradecanoic acid\n'
                 '   - Esters where the chain has unsaturation or other '
                 'modifications\n'
                 '\n'
                 '2. False Negatives: The code is missing legitimate '
                 'tetradecanoate esters, particularly in complex molecules '
                 'where:\n'
                 '   - The chain counting algorithm fails to properly traverse '
                 'the molecule\n'
                 '   - Multiple ester groups are present\n'
                 '   - The tetradecanoate group is part of a larger structure\n'
                 '\n'
                 "Here's how we can improve the code:\n"
                 '\n'
                 '1. Better chain detection:\n'
                 '   - Use SMARTS patterns to specifically match linear '
                 '14-carbon chains\n'
                 '   - Exclude chains with double/triple bonds or '
                 'substitutions\n'
                 '   - Consider both directions from the ester group\n'
                 '\n'
                 '2. Better ester group detection:\n'
                 '   - Specifically look for R-O-C(=O)-R pattern where R is '
                 'not part of a ring\n'
                 '   - Exclude cyclic esters (lactones)\n'
                 '   - Handle cases where the ester oxygen is part of a larger '
                 'structure\n'
                 '\n'
                 '3. Improved validation:\n'
                 '   - Check for presence of modifications on the chain (OH '
                 'groups, double bonds)\n'
                 '   - Verify the chain is terminal (no branches)\n'
                 '   - Handle cases where the ester group is part of a '
                 'phosphate or other complex structure\n'
                 '\n'
                 "Here's my proposed solution:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 68,
    'num_false_positives': 100,
    'num_true_negatives': 27016,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.40476190476190477,
    'recall': 1.0,
    'f1': 0.576271186440678,
    'accuracy': 0.996321365509123,
    'negative_predictive_value': 1.0}