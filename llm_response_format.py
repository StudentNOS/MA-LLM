response_format_schema = {
    "type": "json_schema",
    "json_schema": {
        "name": "pmid_relevance_output_schema",
        "schema": {
            "type": "object",
            "properties": {
                "relevant_pmids": {
                    "description": "List of PMIDs identified as relevant.",
                    "type": "array",
                    "items": {
                        "type": "string",
                        "description": "A single PMID."
                    }
                }
            },
            "required": ["relevant_pmids"],
            "additionalProperties": False
        }
    }
}