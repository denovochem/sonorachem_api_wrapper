# Python wrapper for the SonoraChem API

A python wrapper to access the SonoraChem API.

This API is in pre-alpha. For more details on features, pricing, etc., contact us [here](https://denovochem.com)

## Install

Install directly from this github repo:

```console
pip install git+https://github.com/denovochem/sonorachem_api_wrapper.git
```

## Usage

### Load SonoraChemAPIWrapper and set your API key:
Contact us [here](https://denovochem.com) to obtain an API key
```python
from sonorachem_api import SonoraChemAPIWrapper

api_key = 'API_KEY'
sonorachem_api_wrapper = SonoraChemAPIWrapper(api_key)
```

### Check API usage:

```python
sonorachem_api_wrapper.get_api_usage()
```

### Chemical reaction extraction:

```python
more info coming soon
```

### Forward reaction prediction:

```python
# Single input inference:
reactants = 'C[O-].CC1=C(Cc2c1cccc2)C(O)=O.O.CO.Cl.[Na+].N=C(N)N.CCOC(C)=O'
predictions = sonorachem_api_wrapper.predict_forward_reaction(reactants, sampling_method='top_k', beam_size=5)

# Batch inference:
reactants = [
            'C[O-].CC1=C(Cc2c1cccc2)C(O)=O.O.CO.Cl.[Na+].N=C(N)N.CCOC(C)=O',
            'c1ccncc1.ClCl.ClCCl.O.CS(=O)(=O)CCO.C=C1CC(=O)O1.CC(C)OC(C)C',
            'Cl.OC1CNC1.O=C(Cl)OCc1ccccc1.O=C([O-])[O-].[Na+].[Na+].O.CO.ClC(Cl)Cl.O=S(=O)([O-])[O-].[Na+].[Na+]',
            'CN(C)C=O.Cl.Nc1c(-c2nc3ccc(C4CCNCC4)nc3[nH]2)c(=O)[nH]c2sccc12.CCOC1(O[Si](C)(C)C)CC1.O=C([O-])O.ClCCl.O.[BH3-]C#N.[Cl-].CO.[Na+]',
            'COc1ccccc1C(=O)Cl.Cl[Al](Cl)Cl.O=S(=O)(c1ccccc1)c1ccc(Oc2ccccc2)cc1.ClCCCl.O.[K+].[OH-].O=S(=O)([O-])[O-].[Mg+2]'
            ]
predictions = sonorachem_api_wrapper.batch_predict_forward_reaction(reactants, sampling_method='sampling', temperature=0.3)

print(predictions['output'])
```

### Template-free procedure retrosynthesis prediction:

```python
# Single input inference:
smiles = 'CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5'
predictions = sonorachem_api_wrapper.predict_procedures_retro_template_free(smiles, sampling_method='greedy')

# Batch inference:
smiles = [
          'CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5',
          'CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC',
          'O(CCN(C)C)C(c1ccccc1)c2ccccc2',
          'Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl',
          'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
          ]
predictions = sonorachem_api_wrapper.batch_predict_procedures_retro_template_free(smiles, sampling_method='sampling', temperature=0.7)

print(predictions['output'])
```

### Procedure and condition prediction given reactants and products:

```python
# Single input inference:
reaction_smiles = 'C1N(C)CCN(Cc2c(C(F)(F)F)cc(NC(=O)c3cc(I)c(C)cc3)cc2)C1.c1(C#C)n2ncccc2nc1>>CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5'
predictions = sonorachem_api_wrapper.predict_procedures_given_reactants_products(reaction_smiles, sampling_method='top_k', beam_size=16)

# Batch inference:
reaction_smiles = [
                  'C1N(C)CCN(Cc2c(C(F)(F)F)cc(NC(=O)c3cc(I)c(C)cc3)cc2)C1.c1(C#C)n2ncccc2nc1>>CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5',
                  'C(C(OC1C=C(C(OCC)=O)CC(N=[N+]=[N-])C1NC(=O)C)CC)C>>CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC',
                  'c1cc(C(c2ccccc2)Cl)ccc1.CN(C)CCO>>O(CCN(C)C)C(c1ccccc1)c2ccccc2',
                  'c1(Cl)c(Cl)c(N2CCN(CCCCBr)CC2)ccc1.c12c(cc(O)cc1)NC(=O)CC2>>Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl',
                  'c1(=O)n(C)c2nc[nH]c2c(=O)n1C.CI>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
                  ]
predictions = sonorachem_api_wrapper.batch_predict_procedures_given_reactants_products(reaction_smiles, sampling_method='greedy')

print(predictions['output'])
```

### Purification protocol prediction:

```python
# Single input inference:
reaction_smiles = 'I[Cu].c1cc(P(c2ccccc2)c2ccccc2)ccc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.[Pd].c1c(P(c2ccccc2)c2ccccc2)cccc1.C(N(CC)CC)C.c1c(C)c(C#C)cc(C(Nc2ccc(CN3CCN(C)CC3)c(C(F)(F)F)c2)=O)c1.CCOC(=O)C.O.[Na+].[Cl-]>>CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5'
predictions = sonorachem_api_wrapper.predict_purification_protocols(reaction_smiles, sampling_method='greedy')

# Batch inference:
reaction_smiles = [
                  'I[Cu].c1cc(P(c2ccccc2)c2ccccc2)ccc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.[Pd].c1c(P(c2ccccc2)c2ccccc2)cccc1.C(N(CC)CC)C.c1c(C)c(C#C)cc(C(Nc2ccc(CN3CCN(C)CC3)c(C(F)(F)F)c2)=O)c1.CCOC(=O)C.O.[Na+].[Cl-]>>CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5',
                  'C(C)C(OC1C=C(C(=O)OCC)CC(N=[N+]=[N-])C1NC(=O)C)CC.[H][H].OCC>>CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC',
                  'c1cc(C(c2ccccc2)Cl)ccc1.C(CN(C)C)O.[Na+].[OH-].O.CCOCC>>O(CCN(C)C)C(c1ccccc1)c2ccccc2',
                  'c1(Cl)c(Cl)c(N2CCNCC2)ccc1.C(CCCl)COc1ccc2c(c1)NC(=O)CC2.ClC(Cl)Cl.O.[Na+].[OH-]>>Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl',
                  'c1(=O)n(C)c2nc[nH]c2c(=O)n1C.CI.O=CN(C)C.[K+].[K+].[O-]C([O-])=O.O.ClCCl>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
                  ]
predictions = sonorachem_api_wrapper.batch_predict_purification_protocols(reaction_smiles, sampling_method='sampling', temperature=0.2)

print(predictions['output'])
```

### One-step templated retrosynthesis prediction:

```python
smiles = 'CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5'
predictions = sonorachem_api_wrapper.predict_top_k_retro_templated(smiles, k=16)

print(predictions['output'])
```

### Chemical reaction retrieval:

```python
reaction_smiles = 'I[Cu].c1cc(P(c2ccccc2)c2ccccc2)ccc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.[Pd].c1c(P(c2ccccc2)c2ccccc2)cccc1.C(N(CC)CC)C.c1c(C)c(C#C)cc(C(Nc2ccc(CN3CCN(C)CC3)c(C(F)(F)F)c2)=O)c1.CCOC(=O)C.O.[Na+].[Cl-]>>CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5'

# Search by full reaction
predictions = sonorachem_api_wrapper.retrieve_top_k_similar_reactions_rxn_smiles(reaction_smiles, k=16)

# Search by reactants
predictions = sonorachem_api_wrapper.retrieve_top_k_similar_reactions_reactants(reaction_smiles.split('>>')[0], k=16)

# Search by products
predictions = sonorachem_api_wrapper.retrieve_top_k_similar_reactions_products(reaction_smiles.split('>>')[1], k=16)

print(predictions['output'])
```

## Documentation

Additional documentation coming soon.
