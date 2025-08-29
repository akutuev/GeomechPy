import pytest
from geomechpy.elastic_properties import ElasticProperties, ElasticPropertiesConverter


def test_elastic_properties_converter_with_reference_models() -> None:
    reference_elastic_property_model = ElasticProperties(
        bulk_modulus=40.56,
        youngs_modulus=32.55,
        lames_first_parameter=32.60,
        shear_modulus=11.91,
        poissons_ratio=0.366,
        p_wave_modulus=56.43,
    )

    calculated_elastic_property_results = [
        ElasticPropertiesConverter.from_bulk_and_youngs_modulus(reference_elastic_property_model.bulk_modulus, reference_elastic_property_model.youngs_modulus),
        ElasticPropertiesConverter.from_bulk_and_lames_modulus(reference_elastic_property_model.bulk_modulus, reference_elastic_property_model.lames_first_parameter),
        ElasticPropertiesConverter.from_bulk_and_shear_modulus(reference_elastic_property_model.bulk_modulus, reference_elastic_property_model.shear_modulus),
        ElasticPropertiesConverter.from_bulk_and_poissons_ratio_modulus(reference_elastic_property_model.bulk_modulus, reference_elastic_property_model.poissons_ratio),
        ElasticPropertiesConverter.from_bulk_and_p_wave_modulus(reference_elastic_property_model.bulk_modulus, reference_elastic_property_model.p_wave_modulus),
        ElasticPropertiesConverter.from_youngs_and_lames_modulus(reference_elastic_property_model.youngs_modulus, reference_elastic_property_model.lames_first_parameter),
        ElasticPropertiesConverter.from_youngs_and_shear_modulus(reference_elastic_property_model.youngs_modulus, reference_elastic_property_model.shear_modulus),
        ElasticPropertiesConverter.from_youngs_and_poissons_ratio_modulus(reference_elastic_property_model.youngs_modulus, reference_elastic_property_model.poissons_ratio),
        ElasticPropertiesConverter.from_youngs_and_p_wave_modulus(reference_elastic_property_model.youngs_modulus, reference_elastic_property_model.p_wave_modulus),
        ElasticPropertiesConverter.from_lames_and_shear_modulus(reference_elastic_property_model.lames_first_parameter, reference_elastic_property_model.shear_modulus),
        ElasticPropertiesConverter.from_lames_and_poissons_ratio_modulus(reference_elastic_property_model.lames_first_parameter, reference_elastic_property_model.poissons_ratio),
        ElasticPropertiesConverter.from_lames_and_p_wave_modulus(reference_elastic_property_model.lames_first_parameter, reference_elastic_property_model.p_wave_modulus),
        ElasticPropertiesConverter.from_shear_and_poissons_ratio_modulus(reference_elastic_property_model.shear_modulus, reference_elastic_property_model.poissons_ratio),
        ElasticPropertiesConverter.from_shear_and_p_wave_modulus(reference_elastic_property_model.shear_modulus, reference_elastic_property_model.p_wave_modulus),
        ElasticPropertiesConverter.from_poissons_ratio_and_p_wave_modulus(reference_elastic_property_model.poissons_ratio, reference_elastic_property_model.p_wave_modulus),
    ]

    elastic_properties_test = ElasticPropertiesConverter.from_velocity(3000.0, 2000.0, 2500.0)

    elastic_properties_test.poissons_ratio

    tolerance = 1e-2

    for calculated_elastic_property_result in calculated_elastic_property_results:
        assert reference_elastic_property_model.bulk_modulus == pytest.approx(calculated_elastic_property_result.bulk_modulus, rel=tolerance)
        assert reference_elastic_property_model.youngs_modulus == pytest.approx(calculated_elastic_property_result.youngs_modulus, rel=tolerance)
        assert reference_elastic_property_model.lames_first_parameter == pytest.approx(calculated_elastic_property_result.lames_first_parameter, rel=tolerance)
        assert reference_elastic_property_model.shear_modulus == pytest.approx(calculated_elastic_property_result.shear_modulus, rel=tolerance)
        assert reference_elastic_property_model.poissons_ratio == pytest.approx(calculated_elastic_property_result.poissons_ratio, rel=tolerance)
        assert reference_elastic_property_model.p_wave_modulus == pytest.approx(calculated_elastic_property_result.p_wave_modulus, rel=tolerance)
