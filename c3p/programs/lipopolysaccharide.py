"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:16852 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain:
    - Core oligosaccharide structure
    - Specific sugar units including KDO and heptose
    - Characteristic fatty acid chains
    - Phosphate groups in many cases
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Molecular weight filter
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # LPS components are typically large molecules
        return False, "Molecular weight too low for LPS"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if o_count < 6:  # Need multiple oxygen atoms for glycosidic bonds and hydroxyl groups
        return False, "Too few oxygen atoms for LPS"
    
    if c_count < 15:  # Need sufficient carbons for sugar and lipid components
        return False, "Too few carbon atoms for LPS"

    # Look for sugar patterns more specifically
    pyranose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]([OR0])[CR1]([OR0])[CR1]([OR0])[CR1]")
    furanose_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]([OR0])[CR1]([OR0])[CR1]")
    sugar_matches = len(mol.GetSubstructMatches(pyranose_pattern)) + len(mol.GetSubstructMatches(furanose_pattern))

    # Look for glycosidic linkages more specifically
    glycosidic_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Look for fatty acid chains with hydroxyl groups
    hydroxy_fatty_acid = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H1]")
    fatty_acid_matches = len(mol.GetSubstructMatches(hydroxy_fatty_acid))

    # Look for carboxylic acid groups
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = len(mol.GetSubstructMatches(acid_pattern))

    # Calculate score based on structural features
    score = 0
    
    # Sugar content
    if sugar_matches >= 2:
        score += 2
    elif sugar_matches == 1:
        score += 1

    # Glycosidic linkages
    if glycosidic_matches >= 2:
        score += 2
    elif glycosidic_matches == 1:
        score += 1

    # Fatty acid components
    if fatty_acid_matches >= 1:
        score += 2

    # Acid groups
    if acid_matches >= 1:
        score += 1

    # Phosphate groups
    if p_count >= 1:
        score += 1

    # Ring count
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count >= 3:
        score += 2
    elif ring_count >= 1:
        score += 1

    # Hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches >= 4:
        score += 2
    elif hydroxyl_matches >= 2:
        score += 1

    # Final classification
    if score >= 7:  # Require high confidence for classification
        return True, "Contains characteristic lipopolysaccharide features including sugar units and appropriate linkages"
    else:
        return False, f"Insufficient LPS characteristics (score: {score})"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16412',
                          'name': 'lipopolysaccharide',
                          'definition': 'Liposaccharide natural compounds '
                                        'consisting of a trisaccharide '
                                        'repeating unit (two heptose units and '
                                        'octulosonic acid) with '
                                        'oligosaccharide side chains and '
                                        '3-hydroxytetradecanoic acid units '
                                        '(they are a major constituent of the '
                                        'cell walls of Gram-negative '
                                        'bacteria).',
                          'parents': ['CHEBI:35740', 'CHEBI:65212'],
                          'xrefs': [   'KEGG:C00338',
                                       'PMID:15589368',
                                       'PMID:24186868',
                                       'PMID:24506665',
                                       'PMID:24642373',
                                       'PMID:24659348'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.2192982456140351 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O=C(O)C(O)CCCCCCCCCCCCC(O)CO[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)COC(=O)C '
               'NAME: Ustilagic acid C REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C1OCC2OC(OC(CC(=O)OCC3OC(OC(C1)CCCCCCCCCCCC(OC4OC(C(O)C(C4O)O)CO)CCCCCCCC(O)C)C(O)C(O)C3O)CCCCCCCCCCCCCC(OC5OC(C(O)C(C5O)O)CO)CCCCCCCC(O)C)C(O)C(C2O)O '
               'NAME: Fattiviracin FV-10 REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)C(O)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)C1COC(=O)CCCCCCCCCCCCCCCCC '
               "NAME: alpha,alpha'-Trehalose 4,6-dioleate REASON: CORRECT "
               'Contains characteristic lipopolysaccharide features including '
               'sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O=C1OCC2OC(OC(CC(=O)OCC3OC(OC(C1)CCCCCCCCCCCC(OC4OC(C(O)C(C4O)O)CO)CCCCCC(O)C)C(O)C(O)C3O)CCCCCCCCCCCCCC(OC5OC(C(O)C(C5O)O)CO)CCCCCC(O)C)C(O)C(C2O)O '
               'NAME: Fattiviracin FV-8 REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O1C2(C3(C(C4(C(CC3=O)C(OC4)(C)C)C(O)CC(O)=O)CCC2(C(OC5OC(C(O)C(O)C5O)CO)C=6C=COC6)C)C)C1C(O)=O '
               'NAME: Isolimonic acid glucoside REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O1C2(C3(C(C(C(CC3=O)C(O)(C)C)(C(OC(=O)C)CC(O)=O)C)CCC2(C(OC4OC(C(O)C(O)C4O)CO)C=5C=COC5)C)C)C1C(O)=O '
               'NAME: Nomilinic acid 17-glucoside REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@H](OC(=O)CCC)[C@@H](O[C@@H]([C@H]1OC(=O)C)COC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCCCCCCCCCCCC '
               'NAME: Ustilipid B REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'P(=O)(O[C@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1NC(=O)C[C@H](OC(=O)CCCCCCCCCCC)CCCCCCCCCCC)OC(=O)C[C@H](O)CCCCCCC)CO[C@@H]2O[C@@H]([C@@H](OP(=O)(O)O)[C@@H]([C@H]2NC(=O)CC(=O)CCCCCCCCCCC)OC(=O)C[C@H](O)CCCCCCC)CO)(O)O '
               'NAME: S. salinus M19-40T lipid A REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC[C@H]1O[C@H](O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)COC(=O)[C@@H]([C@H](O)CCCCCCCCCCCCCCC)CCCCCCCCCCCCCC)[C@H](O)[C@H]([C@@H]1O)O)[C@@H]([C@H](O)CCCCCCCCCCCCCCC)CCCCCCCCCCCCCC '
               "NAME: Trehalose-6,6'-dicorynomycolate REASON: CORRECT Contains "
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               '[H][C@]1(O[C@@H](C)[C@H](OC(=O)C(\\C)=C/C)[C@@H](O)[C@H]1OC(=O)[C@@H](C)[C@@H](C)O)O[C@@H]1[C@H]2O[C@]3([H])O[C@H](CO)[C@@H](O)[C@H](OC(=O)[C@@H](C)CC)[C@H]3OC(=O)CCCCCCCCC[C@H](CCCCC)O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O)[C@H]3O[C@]1([H])O[C@H](CO)[C@H]2O '
               'NAME: soldanelline B REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC(CCCCCCCCCCCCCCCC(=O)O[C@@H]1[C@@H](O)[C@H](O[C@@H]([C@H]1O)CO)O[C@H]2O[C@@H]([C@H](O)[C@@H]([C@H]2O)O)CO)CCCCC)C(=O)O '
               'NAME: Emmyguyacin A REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@H](OC(=O)CC)[C@@H](O[C@@H]([C@H]1OC(=O)C)COC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCCCCCCCCCCCC '
               'NAME: Ustilipid E1 REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@@H](O[C@H](COC(=O)C)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C '
               'NAME: Ustilipid E3 REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C1O[C@H]2[C@H](O[C@H]3[C@H](O[C@H](CCCCCCCCCCCCC[C@H]4[C@H]1[C@@H](C(=O)O4)C)C)O[C@H](CO)[C@H]([C@@H]3O)O)O[C@H](CO)[C@H]([C@@H]2O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)O)C)O '
               'NAME: Gobienine C REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C1OCC2OC(OC(CC(=O)OCC3OC(OC(C1)CCCCCCCCCCCCCC(OC4OC(C(O)C(C4OC)O)CO)CCCCCC(=O)C)C(O)C(O)C3O)CCCCCCCCCCCCCCCC(OC5OC(C(O)C(C5OC)O)CO)CCCCCC(OC6OC(C(O)C(C6OC)O)CO)C)C(O)C(C2O)O '
               'NAME: Cycloviracin B2 REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O)[C@@H](NC(=O)CCCCCCCCCCCC[C@H](O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)COC(=O)CCC(=O)O)CC)C '
               'NAME: Ieodoglucomide C REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@H](OC(=O)C)[C@@H](O[C@@H]([C@H]1OC(=O)C)COC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCCCCCCCCCCCC '
               'NAME: Ustilipid E2 REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C1OCC2OC(OC(CC(=O)OCC3OC(OC(C1)CCCCCCCCCCCC(OC4OC(C(O)C(C4O)O)CO)CCCCCC(O)C)C(O)C(O)C3O)CCCCCCCCCCCC(OC5OC(C(O)C(C5O)O)CO)CCCCCC(O)C)C(O)C(C2O)O '
               'NAME: Fattiviracin FV-4 REASON: CORRECT Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O(C1[C@@H](O)[C@H](O)C(O[C@H]1OCCCC)COC(=O)CCCCCCCCCCCCCCC)[C@@H]2OC([C@H](OC(=O)CCC)[C@H](O)C2O)C '
               "NAME: Butyl 4'-O-butanoyl-6-O-hexadecanoyl-neohesperidoside "
               'REASON: CORRECT Contains characteristic lipopolysaccharide '
               'features including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)COC(=O)CCCCCCCCCCCCCC)CCCCCCCCCCCCCC '
               'NAME: '
               '6-O-(6-O-pentadecanoyl-alpha-D-glucopyranosyl)-1-O-pentadecanoyl-alpha-D-glucopyranose '
               'REASON: CORRECT Contains characteristic lipopolysaccharide '
               'features including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)COC(=O)CCCCCCCCCCCCCCC)[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO '
               "NAME: alpha,alpha'-Trehalose 6-palmitate REASON: CORRECT "
               'Contains characteristic lipopolysaccharide features including '
               'sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@@H](O[C@H](COC(=O)C)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCC '
               'NAME: Ustilipid A REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O1C([C@@H](O)C(O)[C@H](O)[C@H]1OCC(OC(=O)CCCCCCCCCCCCCCC)CCCCCCCCCCCCCCCCCC)CO '
               'NAME: '
               '1-O-alpha-D-glucopyranosyl-(2-hexadecanoyloxy)-eicosan-1-ol '
               'REASON: CORRECT Contains characteristic lipopolysaccharide '
               'features including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)C(O)[C@@H](OC1COC(=O)C)OCC(O)CCCCCCCCCCCCC(O)C(O)=O)[C@@H]2OC([C@@H](O)[C@H](O)C2OC(=O)CC(O)CCC)CO '
               'NAME: Ustilagic acid REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C '
               'NAME: Ustilipid D1 REASON: CORRECT Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               'False positives: SMILES: '
               'O1C23OC(CCC2(O)C)CC(CCCC4=NCC(C(CC54C(C(=C(CC5)C6OC(=O)C(=C6)C)C)C=C(C(O)CC7OC1(CC3)CC7)C)C)C)=C '
               'NAME: 13-Demethylspirolide C REASON: WRONGLY CLASSIFIED '
               'Contains characteristic lipopolysaccharide features including '
               'sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O=C1C2(OC(CC(O)CC(O)CC(O)CCCC(O)C=CC(C(C(C(O)C(C(O)C(C(OC3OC(C(OC4OC(C(OC(=O)C(C(OC5OC(C(O)CC5)C)C6=CC=7C(=O)C=C(C)C(C7C=C6)=O)C)CC4)C)C(C3)(O)C)C)CC)C)C)C)OC(CCC(C(CC(CCCC(CCC(C1O)C)O)O)O)C)=O)C)CC(C2)O)O '
               'NAME: Iseolide B REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@H]([C@H](O)C[C@H](O)[C@@H]([C@@H](OC)C[C@H]2O[C@@H](C=CC2)C[C@@H](O)CC=C(C(=O)O[C@@H]([C@H]([C@@H](C[C@@H]([C@@H]([C@H](C[C@H]3O[C@H](C[C@H](CC=C1C)O)C=CC3)OC)C)O)O)C)[C@H]([C@@H](O)[C@H](CC[C@@H]4O[C@H](C[C@H](C4)OC)C)C)C)C)C)C)[C@H]([C@@H](O)[C@H](CC[C@@H]5O[C@H](C[C@H](C5)OC)C)C)C '
               'NAME: Misakinolide A REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC[C@@H](O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: '
               '1-(3-O-sulfo-beta-D-galactosyl)-N-[(2R)-2-hydroxystearoyl]sphingosine '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'OC12C3C(C4(C(CC3)CC(OC5OC(C(O)C(O)C5OC6OC(C(O)C(O)C6O)CO)CO)CC4)C)CCC1(C(CC2)C=7COC(=O)C7)C '
               'NAME: '
               '3-[3-[4,5-Dihydroxy-6-(hydroxymethyl)-3-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-14-hydroxy-10,13-dimethyl-1,2,3,4,5,6,7,8,9,11,12,15,16,17-tetradecahydrocyclopenta[a]phenanthren-17-yl]-2H-furan-5-one '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(O[C@H]1[C@H](O[C@H](COC(=O)C)[C@H]([C@@H]1O)O)O[C@@H]2O[C@H]([C@H](O)[C@H]([C@@H]2O)O)CO)CCCCCCCCCCCCCCCCC '
               'NAME: 6-0-acetyl2-0-octadecanoy 1-alpha ,alpha-D-trehalose '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](OCCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O '
               'NAME: bhos#32-CoA REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@@]2(C[C@@H](C1)O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C([O-])=O)[H])[H])(CC[C@]4([H])[C@@H](O)C)[H])C)[H])C '
               'NAME: 5beta-pregnane 3alpha,20alpha-diol '
               '3-O-(beta-D-glucuronide)(1-) REASON: WRONGLY CLASSIFIED '
               'Contains characteristic lipopolysaccharide features including '
               'sugar units and appropriate linkages\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxydocosanoyl)-1-O-beta-D-glucosyl-15-methylhexadecasphing-4-enine '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'S(CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]1O[C@@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H]([C@@H]1OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O)C(\\C=C\\CCCCCCCCCCCCCCCO[C@H]4[C@@H](C[C@H]([C@@H](O4)C)O)O)=O '
               'NAME: oscr#31-CoA REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O1C(OC2C(O)C(O)C(OC2CO)OCC(NC(=O)CCCCCCCCCCCCCCC)C(O)/C=C/CCCCCCCCCCCCC)C(O)C(O)C(OC3OC(C(O)C(O)C3O)CO)C1CO '
               'NAME: Trihexosylceramide (d18:1/16:0) REASON: WRONGLY '
               'CLASSIFIED Contains characteristic lipopolysaccharide features '
               'including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'C1[C@](OC(C[C@@H](C(C([C@@H]([C@H]([C@H](CCC[C@@]2([C@]1(O2)[H])C)C)O)C)=O)(C)C)O)=O)(/C(=C/C=3N=C(SC3)C)/C)[H] '
               'NAME: epothilone B REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: '
               "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/24:0) "
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O[C@@]12[C@]3([C@@]([C@@]4([C@](O)(CC3)C[C@@H](O[C@@H]5O[C@H]([C@H](O)[C@@H](O)[C@H]5O)C)C[C@H]4O)CO)([C@H](O)C[C@@]1([C@H](CC2)C=6COC(=O)C6)C)[H])[H].O '
               'NAME: 11,15>]Heptadec-14-yl]-5-hydrofuran-2-one, hydrate '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               '[H][C@@]12CC=C3C[C@H](CC[C@]3(C)[C@@]1([H])CC[C@]1(C)C(=O)CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O '
               'NAME: dehydroisoandrosterone 3-glucuronide REASON: WRONGLY '
               'CLASSIFIED Contains characteristic lipopolysaccharide features '
               'including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'C1[C@]2([C@]3([C@@](C4=C(C=C(O)C=C4)CC3)(CC[C@@]2(C(=O)[C@@H]1O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C(O)=O)C)[H])[H])[H] '
               'NAME: 16alpha-hydroxyestrone 16-O-(beta-D-glucuronide) REASON: '
               'WRONGLY CLASSIFIED Contains characteristic lipopolysaccharide '
               'features including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'O1C(O[C@@H]2[C@H](O)[C@H](O[C@@H]([C@H]2O)CO)OCC(O)COCCCCCCCCCCCC(CC)C)[C@H](O)[C@@H](O)[C@@H]([C@H]1CO)O '
               'NAME: GGL 2 REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O1C2C(C3(C(C4C(CC3)C5(C(=CC4)CC(OC6OC(C(OC7OC(C(O)C(OC8OC(C(O)C(O)C8O)CO)C7O)CO)C(O)C6OC9OC(C(O)C(O)C9O)C)CO)CC5)C)C2)C)C(C1(O)CCC(COC%10OC(C(O)C(O)C%10O)CO)C)C '
               'NAME: Trigoneoside XIIIa REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O=C1[C@@H](OC)[C@@](O)([C@@](O)(COC(=O)[C@H](O)C[C@@]23O[C@@](C=CC=C(C)C)(C)[C@@H](C2)C(CC3)=C)CC1)[C@]4(O[C@H]4CC=C(C)C)C '
               'NAME: Pseudallicin A REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](OS([O-])(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-docosanoyl-(3-O-sulfo-D-galactosyl)sphingosine(1-) '
               'REASON: WRONGLY CLASSIFIED Contains characteristic '
               'lipopolysaccharide features including sugar units and '
               'appropriate linkages\n'
               ' * SMILES: '
               'C[C@H](CCCCCCCCCCCCCC\\C=C\\C(O)=O)O[C@@H]1O[C@@H](C)[C@H](O)C[C@H]1O '
               'NAME: ascr#33 REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)[C@@H]1OCC(NC(=O)CCCCCCCCCCCCCCCCCCC)C(O)/C=C\\CCCCCCCCCCCCC)CO '
               'NAME: Galactosylceramide (d18:1/20:0) REASON: WRONGLY '
               'CLASSIFIED Contains characteristic lipopolysaccharide features '
               'including sugar units and appropriate linkages\n'
               ' * SMILES: '
               'CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(O)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)[NH+](C)C)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@H]1C '
               'NAME: erythromycin D(1+) REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(C(O)CC(C5)(C)C)C)=CC3)C)(CC2)C)(CC1)C)(C)C)C6OC(C(O)C(O)C6OC7OC(C(O)C(O)C7O)CO)C(O)=O '
               'NAME: Azukisaponin I REASON: WRONGLY CLASSIFIED Contains '
               'characteristic lipopolysaccharide features including sugar '
               'units and appropriate linkages\n'
               ' * SMILES: '
               'CC(=O)OC(C)(C)[C@H]1CC[C@](O1)(C)[C@H](CC[C@@](C)([C@]2(CC[C@@]3([C@@](O2)(CC[C@@](O3)([C@@]4(CC[C@H](C(O4)(C)C)Br)C)[H])C)[H])[H])O)O '
               'NAME: thyrsiferyl 23-acetate REASON: WRONGLY CLASSIFIED '
               'Contains characteristic lipopolysaccharide features including '
               'sugar units and appropriate linkages\n'
               'False negatives: SMILES: '
               'O1[C@H](O[C@H]([C@H](O)[C@@H](O)C=O)[C@H](O)C(O)=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1C(O)=O '
               'NAME: 4-O-alpha-Galactopyranuronosyl-galacturonic acid REASON: '
               'MISSED Insufficient characteristic LPS features\n'
               ' * SMILES: '
               'CCCCCCCCCC(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H](C[C@H](OC(C)=O)C(\\C)=C\\CO[C@H]([C@H](O)CO)[C@H](O)[C@H](O)CO)[C@](C)(O)CCC=C(C)C '
               'NAME: Cosmosporaside D, (rel)- REASON: MISSED Insufficient '
               'characteristic LPS features\n'
               ' * SMILES: '
               'CC(C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](OC(=O)CC(C)(O)CC(=O)O[C@H](COc2c3ccoc3cc3oc(=O)ccc23)C(C)(C)O)[C@H]1O '
               'NAME: citrusoside D REASON: MISSED Insufficient characteristic '
               'LPS features\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@H](O)[C@H](O)[C@@H]1O)C(O)=O)[C@H]([C@@H](O)CCCCCCCC(O)=O)C/C=C\\CCCCC '
               'NAME: 9-Hydroxy-10-O-D-glucuronoside-12Z-octadecenoate REASON: '
               'MISSED Insufficient characteristic LPS features\n'
               ' * SMILES: '
               'O1[C@](O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)([C@@]3(OC(=O)C4C(C(C=5C4=CC(O)=C(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)COC(=O)CC7C=8C(NC7=O)=CC=CC8)C5)C9=CC(OC)=C(O)C=C9)C(OC[C@]1([C@@H]3O)[H])=O)[H])CO '
               'NAME: Secaloside A REASON: MISSED Insufficient characteristic '
               'LPS features\n'
               ' * SMILES: '
               'CC(C)O[C@@H]1O[C@H](COC(=O)CC(C)(O)CC(=O)O[C@H](CC\\C(C)=C\\COc2c3ccoc3cc3oc(=O)ccc23)C(C)(C)O)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: citrusoside B REASON: MISSED Insufficient characteristic '
               'LPS features\n'
               ' * SMILES: '
               'CC(C)O[C@@H]1O[C@H](COC(=O)CC(C)(O)CC(=O)O[C@H](COc2c3ccoc3cc3oc(=O)ccc23)C(C)(C)O)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: citrusoside C REASON: MISSED Insufficient characteristic '
               'LPS features\n'
               ' * SMILES: '
               'O1C(OC2=C(O)C=C(CC(O)CCC(O)=O)C=C2)C(O)C(O)C(O)C1C(O)=O NAME: '
               '6-[4-(4-carboxy-2-hydroxybutyl)-2-hydroxyphenoxy]-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: MISSED Insufficient characteristic LPS features\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@H](OC(=O)C)[C@@H](O[C@@H]([C@H]1OC(=O)C)CO)OC[C@@H](O)[C@@H](O)CO)CCCCCCCCCCCCC '
               'NAME: Ustilipid F3 REASON: MISSED Insufficient characteristic '
               'LPS features\n'
               ' * SMILES: O1C(C(O)C(OC(=O)CCCCC)C(OC(=O)CCCCC)C1O)CO NAME: '
               '2,3-di-0-hexanoyl-alpha-glucopyranose REASON: MISSED '
               'Insufficient characteristic LPS features\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'Insufficient LPS '
                                               'characteristics (score: 5)'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'Molecular weight too low for '
                                               'LPS'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'Too few oxygen atoms for LPS'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'Too few oxygen atoms for LPS'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'Molecular weight too low for '
                                               'LPS'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'Too few oxygen atoms for LPS'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'Too few oxygen atoms for LPS'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'Insufficient LPS '
                                               'characteristics (score: 3)'},
                                 {   'smiles': 'OC([C@@H](O)CC1=CC=2NC=C(C2C=C1)CC#N)(C)C',
                                     'name': 'Streptoprenylindole A',
                                     'reason': 'Molecular weight too low for '
                                               'LPS'},
                                 {   'smiles': 'C[C@H]1[C@H]2[C@H](C[C@@H]3[C@@]2(CC[C@H]4[C@H]3CC[C@@H]5[C@@]4(CC[C@@H](C5)O)C)C)O[C@@]16CCC(CN6)C',
                                     'name': 'LSM-6471',
                                     'reason': 'Too few oxygen atoms for LPS'}],
    'sample_false_negatives': [   {   'smiles': 'O1[C@H](O[C@H]([C@H](O)[C@@H](O)C=O)[C@H](O)C(O)=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1C(O)=O',
                                      'name': '4-O-alpha-Galactopyranuronosyl-galacturonic '
                                              'acid',
                                      'reason': 'Too few carbon atoms for LPS'},
                                  {   'smiles': 'CCCCCCCCCC(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H](C[C@H](OC(C)=O)C(\\C)=C\\CO[C@H]([C@H](O)CO)[C@H](O)[C@H](O)CO)[C@](C)(O)CCC=C(C)C',
                                      'name': 'Cosmosporaside D, (rel)-',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 6)'},
                                  {   'smiles': 'O=C(O[C@@H]1[C@H](OC(=O)C)[C@@H](O[C@@H]([C@H]1OC(=O)C)CO)OC[C@@H](O)[C@@H](O)CO)CCCCCCCCCCCCC',
                                      'name': 'Ustilipid F3',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 6)'},
                                  {   'smiles': 'O1C(C(O)C(OC(=O)CCCCC)C(OC(=O)CCCCC)C1O)CO',
                                      'name': '2,3-di-0-hexanoyl-alpha-glucopyranose',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 5)'},
                                  {   'smiles': 'O1C([C@@H](O)C(O)C(O)[C@@H]1O)COC(=O)C([C@@H](O)CO)=C',
                                      'name': '6-O-(3R,4-dihydroxy-2-methylene-butanoyl)-beta-D-glucopyranose',
                                      'reason': 'Molecular weight too low for '
                                                'LPS'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OCCC(C(O)=O)=C)CO',
                                      'name': '2-methylidene-4-[(2r,3r,4s,5s,6r)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxybutanoic '
                                              'acid',
                                      'reason': 'Molecular weight too low for '
                                                'LPS'},
                                  {   'smiles': 'O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC(C(O)(CCC=C(C)C)C)CC(OC(=O)C)C(=CCO[C@@H]([C@H](O)[C@H](O)CO)[C@H](O)CO)C)CCCCCCCCC',
                                      'name': 'Cosmosporaside D',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 6)'},
                                  {   'smiles': 'CCOC(=O)[C@H]1CC=CC[C@@H](C(=O)O[C@@H]2[C@H]([C@@H]([C@H](O[C@@H]2OC)CO)OC(=O)NC3=CC=C(C=C3)Br)OC1=O)C(=O)OCC',
                                      'name': 'LSM-36740',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 3)'},
                                  {   'smiles': 'CCOC(=O)[C@H]1CC=CC[C@@H](C(=O)O[C@@H]2[C@H]([C@@H]([C@H](O[C@@H]2OC)CO)OC(=O)NC3=CC=C(C=C3)OC)OC1=O)C(=O)OCC',
                                      'name': 'LSM-36728',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 3)'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC(=O)/C=C/C2=CC=CC=C2)COC(=O)C(C(O)CO)=C',
                                      'name': '[(2R,3S,4S,5R,6S)-3,4,5-trihydroxy-6-[(E)-3-phenylprop-2-enoyl]oxyoxan-2-yl]methyl '
                                              '3,4-dihydroxy-2-methylidenebutanoate',
                                      'reason': 'Insufficient LPS '
                                                'characteristics (score: 6)'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. False Positives Analysis:\n'
                 '- Many complex glycosides and glycoconjugates were '
                 'incorrectly classified as LPS\n'
                 '- The program was too lenient in its classification, mainly '
                 'looking for sugar units and fatty acid chains\n'
                 '- Many steroid glycosides and other glycoconjugates have '
                 "similar patterns but aren't LPS\n"
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- Several legitimate LPS components were missed\n'
                 '- The criteria for KDO (keto-deoxy-octulosonic acid) was too '
                 'strict\n'
                 '- The requirement for heptose units was too rigid\n'
                 '\n'
                 '3. Key Issues:\n'
                 '- The structural requirements were not specific enough\n'
                 '- The scoring system (`lps_features`) was too simple\n'
                 '- The fatty acid pattern matching was too general\n'
                 '\n'
                 '4. Proposed Improvements:\n'
                 '- Make pattern matching more specific to LPS core '
                 'structures\n'
                 '- Add better distinction between LPS and other '
                 'glycoconjugates\n'
                 '- Implement more precise sugar unit identification\n'
                 '- Add molecular weight constraints\n'
                 '- Consider the presence of phosphate groups\n'
                 '- Look for specific linkage patterns characteristic of LPS\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 71,
    'num_false_positives': 100,
    'num_true_negatives': 823,
    'num_false_negatives': 32,
    'num_negatives': None,
    'precision': 0.4152046783625731,
    'recall': 0.6893203883495146,
    'f1': 0.5182481751824818,
    'accuracy': 0.8713450292397661,
    'negative_predictive_value': 0.9625730994152046}