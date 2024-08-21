def round_dict_content(
    input: dict[str, float], precision: int, plus: int = 0
) -> dict[str, float]:
    return {k: round(v, precision + plus) for k, v in input.items()}
