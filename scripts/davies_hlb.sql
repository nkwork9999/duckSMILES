-- Load this reference in a separate DuckDB connection with .read or execute().
CREATE OR REPLACE TEMP MACRO davies_hlb_ref(counts) AS (
    WITH groups(group_id, signed_milli) AS (
        VALUES ('sodium_sulfate', 38700),
               ('potassium_carboxylate', 21100),
               ('sodium_carboxylate', 19100),
               ('tertiary_amine', 9400),
               ('sorbitan_ester', 6800),
               ('free_ester', 2400),
               ('carboxylic_acid', 2100),
               ('free_hydroxyl', 1900),
               ('ether', 1300),
               ('sorbitan_hydroxyl', 500),
               ('hydrocarbon', -475),
               ('ethylene_oxide', 330),
               ('propylene_oxide', -150)
    )
    SELECT CASE
        WHEN count(*) FILTER (WHERE signed_milli IS NULL OR entry.value IS NULL
                              OR entry.value < 0) = 0
             AND count(*) FILTER (WHERE entry.value > 0) > 0
        THEN CAST(7000 + sum(entry.value::HUGEINT * signed_milli) AS DOUBLE) / 1000.0
    END
    FROM unnest(map_entries(counts::MAP(VARCHAR, BIGINT))) AS m(entry)
    LEFT JOIN groups ON encode(entry.key) = encode(groups.group_id)
);
