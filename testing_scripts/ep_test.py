from geomechpy.elastic_properties import ElasticProperties, ElasticPropertiesConverter


result1 = ElasticPropertiesConverter.from_bulk_and_youngs_modulus(19.7184, 33.2748)
result2 = ElasticPropertiesConverter.from_bulk_and_lame_modulus(19.7184, 10.61851)
