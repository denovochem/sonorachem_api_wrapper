# Python wrapper for the SonoraChem API

A python wrapper to access the SonoraChem API.

This API is in pre-alpha. For more details on features, pricing, etc., contact us [here](https://denovochem.com)

## Install

Install directly from this github repo:

```console
pip install git+https://github.com/denovochem/sonorachem_api_wrapper.git
```

## Usage

### Load SonoraChemAPIWrapper and set an API key:

```python
api_key = 'API_KEY'
from sonorachem_api import SonoraChemAPIWrapper
sonorachem_api_wrapper = SonoraChemAPIWrapper(api_key)
```

### Check API usage:

```python
sonorachem_api_wrapper.get_api_usage()
```

### Forward reaction prediction:

```python
reactants = ['COc1ccccc1C(=O)Cl.Cl[Al](Cl)Cl.O=S(=O)(c1ccccc1)c1ccc(Oc2ccccc2)cc1.ClCCCl.O.[K+].[OH-].O=S(=O)([O-])[O-].[Mg+2]']
predictions = sonorachem_api_wrapper.predict_forward_reaction(reactants)
```

### Template-free procedure retrosynthesis prediction:

```python
smiles = ['CC1=C(C=C(C=C1)C(=O)NC2=CC(=C(C=C2)CN3CCN(CC3)C)C(F)(F)F)C#CC4=CN=C5N4N=CC=C5']
predictions = sonorachem_api_wrapper.predict_procedures_retro_template_free(smiles)
```

### Procedure and condition prediction given reactants and products:

```python
reaction_smiles = ['c1cc(C(c2ccccc2)Cl)ccc1.CN(C)CCO>>O(CCN(C)C)C(c1ccccc1)c2ccccc2']
predictions = sonorachem_api_wrapper.predict_procedures_given_reactants_products(reaction_smiles)
```

### Purification protocol prediction:

```python
reaction_smiles = ['c1(Cl)c(Cl)c(N2CCNCC2)ccc1.C(CCCl)COc1ccc2c(c1)NC(=O)CC2.ClC(Cl)Cl.O.[Na+].[OH-]>>Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl']
predictions = sonorachem_api_wrapper.predict_purification_protocols(reaction_smiles)
```

## Documentation

Additional documentation coming soon.

