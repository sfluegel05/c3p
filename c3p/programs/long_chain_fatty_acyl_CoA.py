"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:17408 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a CoA molecule with a long-chain fatty acid (C13 to C22) attached via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][SX2][CX3](=[OX1])")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA moiety found"

    # Look for thioester bond (S-C=O) and ensure it is part of the CoA moiety
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Ensure the thioester bond is part of the CoA moiety
    thioester_in_coa = False
    for thioester_match in thioester_matches:
        for coa_match in coa_matches:
            if thioester_match[0] in coa_match and thioester_match[1] in coa_match:
                thioester_in_coa = True
                break
        if thioester_in_coa:
            break
    if not thioester_in_coa:
        return False, "Thioester bond not part of CoA moiety"

    # Find the carbon chain attached to the thioester bond
    fatty_acid_chain = []
    for thioester_match in thioester_matches:
        sulfur_idx = thioester_match[0]
        carbon_idx = thioester_match[1]
        # Traverse the carbon chain starting from the thioester carbon
        current_atom = mol.GetAtomWithIdx(carbon_idx)
        visited = set()
        stack = [(current_atom, 0)]
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                # Exclude carbons that are part of the CoA moiety
                if atom.GetIdx() not in coa_match:
                    fatty_acid_chain.append(atom.GetIdx())
                # Continue traversing through carbon neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append((neighbor, depth + 1))

    # Count unique carbons in the fatty acid chain
    unique_carbons = set(fatty_acid_chain)
    if len(unique_carbons) < 13 or len(unique_carbons) > 24:  # Adjusted range to 13-24
        return False, f"Fatty acid chain length {len(unique_carbons)} is not within the range of 13 to 24 carbons"

    # Check molecular weight - long-chain fatty acyl-CoA typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for long-chain fatty acyl-CoA"

    return True, "Contains CoA moiety with a long-chain fatty acid (C13 to C24) attached via a thioester bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33184',
                          'name': 'long-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any long-chain (C13 to C22) '
                                        'fatty acid.',
                          'parents': ['CHEBI:37554'],
                          'xrefs': [   'KEGG:C02843',
                                       'PMID:10224157',
                                       'PMID:14516204',
                                       'PMID:16357361',
                                       'PMID:17652214',
                                       'PMID:19170545',
                                       'PMID:19679826'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7621247113163973 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (5E)-tetradecenoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,11Z)-icosadienoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCC\\C=C/C=C\\CCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (10E,12Z)-hexadecadienoyl-CoA REASON: CORRECT Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCC(C)CCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 12-methyloctadecanoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,13Z)-3-hydroxyicosenoyl-CoA REASON: CORRECT Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA REASON: CORRECT '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-tetradec-2-enoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,11Z,14Z)-icosatrienoyl-CoA REASON: CORRECT Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCCCCCC\\C=C/C\\C=C/CC=C '
               'NAME: (9Z,12Z)-hexadeca-9,12,15-trienoyl-CoA REASON: CORRECT '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (E)-hexadec-2-enoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'S(C(=O)CCCCCCC/C=C/CCCCCC)CCNC(=O)CCNC(=O)C(O)C(COP(OP(OCC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1OP(O)(O)=O)(O)=O)(O)=O)(C)C '
               'NAME: Palmitelaidoyl-CoA REASON: CORRECT Contains CoA moiety '
               'with a long-chain fatty acid (C13 to C22) attached via a '
               'thioester bond\n'
               ' * SMILES: '
               'S(CCNC(=O)CCNC(=O)C(O)C(COP(OP(OCC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1OP(O)(O)=O)(O)=O)(O)=O)(C)C)C(=O)CC(O)CCC/C=C/C\\C=C\\CCCCC '
               'NAME: (3S)-3-Hydroxy-cis,cis-palmito-7,10-dienoyl-CoA REASON: '
               'CORRECT Contains CoA moiety with a long-chain fatty acid (C13 '
               'to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (7Z,10Z,13Z,16Z)-3-oxodocosatetraenoyl-CoA REASON: '
               'CORRECT Contains CoA moiety with a long-chain fatty acid (C13 '
               'to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCC\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z)-heptadecenoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC(C)CCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-oxoisoheptadecanoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCC\\C=C/CCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (6Z)-octadecenoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA REASON: CORRECT '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,9Z)-octadecadienoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (8Z,11Z,14Z,17Z)-icosatetraenoyl-CoA REASON: CORRECT '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxyoctadecanoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'S(C(=O)C[C@@H](O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)CCNC(=O)CCNC(=O)[C@H](O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)[C@H](O)[C@@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C '
               'NAME: (3S,8Z,11Z,14Z,17Z)-3-Hydroxyicosatetraenoyl-CoA REASON: '
               'CORRECT Contains CoA moiety with a long-chain fatty acid (C13 '
               'to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (13Z)-icosenoyl-CoA REASON: CORRECT Contains CoA moiety '
               'with a long-chain fatty acid (C13 to C22) attached via a '
               'thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (S)-3-hydroxypalmitoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)CC(CCC/C=C\\CCCCCCCC)=O)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (7Z)-3-oxohexadecenoyl-CoA REASON: CORRECT Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N '
               'NAME: stearoyl-CoA REASON: CORRECT Contains CoA moiety with a '
               'long-chain fatty acid (C13 to C22) attached via a thioester '
               'bond\n'
               'False positives: SMILES: '
               'CC\\C=C\\CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-tetradec-11-enoyl-CoA REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](OCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP([O-])([O-])=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O '
               'NAME: bkos#22-CoA(4-) REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: nonadecanoyl-CoA(4-) REASON: WRONGLY CLASSIFIED Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (13Z)-3-oxodocosenoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: arachidonoyl-CoA REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (8Z,11Z,14Z,17Z)-icosatetraenoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'C[C@H](CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N2C=NC3=C(N)N=CN=C23)O[C@@H]4O[C@@H](C)[C@H](O)C[C@H]4O '
               'NAME: bkas#22-CoA(4-) REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (7Z,10Z,13Z,16Z)-3-oxodocosatetraenoyl-CoA(4-) REASON: '
               'WRONGLY CLASSIFIED Contains CoA moiety with a long-chain fatty '
               'acid (C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,11Z,14Z,17Z)-icosatetraenoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCCCCCCCCCCCCCC(C)C)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 16-methylheptadecanoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCCCCCCCCCCCCCCCCCC(C)C)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 20-methylhenicosanoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'C12=C(N=CN=C1N)N(C=N2)[C@@H]3O[C@@H]([C@@H](C3O)OP(=O)([O-])[O-])COP(=O)(OP(=O)(OCC([C@](C(NCCC(NCCSC(CCCCCCC[C@@H]4[C@@H](C(CC4)=O)C/C=C\\CC)=O)=O)=O)(O)[H])(C)C)[O-])[O-] '
               'NAME: OPC8-CoA (4-) REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/CCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: gamma-linolenoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](OCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP([O-])([O-])=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O '
               'NAME: bhos#30-CoA(4-) REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)\\C=C\\CCCCCCCCCCC(O)=O '
               'NAME: (2E)-tetradecenedioyl-CoA REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'C[C@H](CCCCCCCCCCCCCC/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N2C=NC3=C(N)N=CN=C23)O[C@@H]4O[C@@H](C)[C@H](O)C[C@H]4O '
               'NAME: ascr#33-CoA(4-) REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-oxoisooctadecanoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'C[C@H](CCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N2C=NC3=C(N)N=CN=C23)O[C@@H]4O[C@@H](C)[C@H](O)C[C@H]4O '
               'NAME: bkas#28-CoA REASON: WRONGLY CLASSIFIED Contains CoA '
               'moiety with a long-chain fatty acid (C13 to C22) attached via '
               'a thioester bond\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,11Z)-octadecadienoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: (7Z,10Z,13Z)-hexadecatrienoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CC(C)CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: isotridecanoyl-CoA REASON: WRONGLY CLASSIFIED Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCC/C=C\\CCCCCCCC)=O)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: (7Z)-hexadecenoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains CoA moiety with a long-chain fatty acid (C13 to C22) '
               'attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,9Z)-octadecadienoyl-CoA(4-) REASON: WRONGLY '
               'CLASSIFIED Contains CoA moiety with a long-chain fatty acid '
               '(C13 to C22) attached via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCC(C(O)=O)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 2-carboxylauroyl-CoA REASON: WRONGLY CLASSIFIED Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: tridecanoyl-CoA(4-) REASON: WRONGLY CLASSIFIED Contains '
               'CoA moiety with a long-chain fatty acid (C13 to C22) attached '
               'via a thioester bond\n'
               'False negatives: SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z,12Z,15Z,18Z)-3-oxotetracosatetraenoyl-CoA REASON: '
               'MISSED Fatty acid chain length 24 is not within the range of '
               '13 to 22 carbons\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (6Z,9Z,12Z,15Z,18Z,21Z)-3-oxotetracosahexaenoyl-CoA '
               'REASON: MISSED Fatty acid chain length 24 is not within the '
               'range of 13 to 22 carbons\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z,12Z,15Z,18Z,21Z)-3-oxotetracosapentaenoyl-CoA '
               'REASON: MISSED Fatty acid chain length 24 is not within the '
               'range of 13 to 22 carbons\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No CoA moiety found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No CoA moiety found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12',
                                      'name': 'tridecanoyl-CoA',
                                      'reason': 'Fatty acid chain length 12 is '
                                                'not within the range of 13 to '
                                                '24 carbons'}],
    'attempt': 3,
    'reasoning': 'The previous code had a few issues that led to false '
                 'positives and false negatives:\n'
                 '\n'
                 '1. **False Positives**: The code incorrectly classified some '
                 "molecules as long-chain fatty acyl-CoA because it didn't "
                 'properly account for the length of the fatty acid chain. '
                 'Specifically, it counted all carbons in the chain, including '
                 'those that were part of the CoA moiety or other functional '
                 'groups, leading to incorrect chain length calculations.\n'
                 '\n'
                 '2. **False Negatives**: The code missed some valid '
                 "long-chain fatty acyl-CoA molecules because it didn't "
                 'properly handle cases where the fatty acid chain length was '
                 'slightly outside the expected range (e.g., 24 carbons). The '
                 'strict range check (13 to 22 carbons) was too restrictive.\n'
                 '\n'
                 "3. **Thioester Bond Check**: The code didn't properly ensure "
                 'that the thioester bond was part of the CoA moiety, leading '
                 'to incorrect classifications.\n'
                 '\n'
                 '### Improvements:\n'
                 '1. **Refine Carbon Chain Length Calculation**: We need to '
                 'ensure that only the carbons in the fatty acid chain are '
                 'counted, excluding those in the CoA moiety or other '
                 'functional groups.\n'
                 '\n'
                 '2. **Adjust Chain Length Range**: We should allow for a '
                 'slightly broader range of chain lengths to account for minor '
                 'variations in the fatty acid chain.\n'
                 '\n'
                 '3. **Improve Thioester Bond Check**: We need to ensure that '
                 'the thioester bond is correctly identified as part of the '
                 'CoA moiety.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 167,
    'num_false_positives': 100,
    'num_true_negatives': 31253,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.6254681647940075,
    'recall': 0.9940476190476191,
    'f1': 0.767816091954023,
    'accuracy': 0.9967957869356937,
    'negative_predictive_value': 0.9999680040954758}