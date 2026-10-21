#!/usr/bin/env python3
"""Native SQL checks for 16 genomic, IRT, SEM and degradation operations."""
import argparse
import json
import math
from pathlib import Path

from check_jmp_continuation import Runner, equal, expression, field, sql_list, strict_json

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb",required=True)
    parser.add_argument("--extension",required=True,type=Path)
    args=parser.parse_args()
    runner=Runner(args.duckdb,args.extension)
    apis={g:json.loads((ROOT/f"tasks/jmp-{g}-api.json").read_text()) for g in ["kinship","irt","sem","degradation"]}
    operations=[o for api in apis.values() for o in api["operations"]]
    assert len(operations)==len({o["operation"]for o in operations})==16
    names=",".join("'"+o["sql_name"]+"'" for o in operations)
    catalog=runner.query(f"SELECT function_name,description,examples,parameter_types FROM duckdb_functions() WHERE function_name IN({names});")
    assert len(catalog)==16 and all(r["description"]and r["examples"]and r["parameter_types"]==["DOUBLE[][]","DOUBLE[]"]for r in catalog)
    examples=[o["sql_example"]for o in operations]
    numeric_fields=0
    for o,result in zip(operations,runner.many(examples)):
        assert o["known_value"],o["operation"]
        for key,wanted in o["known_value"].items():
            equal(field(result,key),wanted,o["operation"]+"."+key,5e-5)
            numeric_fields+=1
    reuse=[]
    def saved(fit,predict,queries,expected,key="predict_parameters",prefix=()):
        fit_call=next(o["sql_example"]for o in operations if o["operation"]==fit)
        params=f"json_extract(value,'$.{key}')::DOUBLE[]"
        if prefix:
            params=f"list_concat({sql_list(prefix)},{params})"
        sql=f"(WITH model AS(SELECT {fit_call} AS value) SELECT chem_{predict}({sql_list(queries)},{params}) FROM model)"
        reuse.append((predict,sql,expected))
    saved("kinship_model","kinship_predict",[[1,1,0,0,0,0],[0]*6],{"predictions":[1.25,5]})
    saved("irt_fit","irt_predict",[[0]],{"probabilities":[[.5,.5]]},"model_parameters")
    saved("irt_fit","irt_score",[[0,0],[1,1]],{"posterior_means":[-.7042799208115517,.7042799208115517]},"model_parameters")
    means=[row[-1]for row in next(o["fixture"]["data"]for o in operations if o["operation"]=="cfa_model")]
    saved("cfa_model","cfa_scores",[means],{"scores":[[0,0]]},"score_parameters",[0])
    saved("degradation_path","degradation_predict",[[2,1,6]],{"fitted":[6],"threshold_times":[2]},prefix=[.95,0])
    saved("destructive_degradation","destructive_predict",[[2,0,6]],{"medians":[6],"means":[6],"lower_failure_probabilities":[.5]})
    saved("cumulative_damage_fit","cumulative_damage_predict",[[2,1,.5]],{"effective_ages":[1],"survivals":[math.exp(-1)],"hazards":[.5],"quantiles":[2*math.log(2)]})
    for predict,sql,expected in reuse:
        result=runner.one(sql)
        for key,wanted in expected.items():
            equal(field(result,key),wanted,"saved "+predict+"."+key,5e-5)
            numeric_fields+=1
    ram=runner.one(expression("sem_ram",[[0,2,1,0,0],[2,0,0,1,0]]))
    equal(ram["covariance"],[[5/9,4/9],[4/9,5/9]],"invertible RAM covariance")
    assert not ram["dynamically_stable"],ram
    invalid=[
        expression("kinship_predict",[[1]],[]),expression("kinship_model",[[1,1,2],[2,2,1]],[0]),
        expression("irt_predict",[[0]],[1,1,1,0,41,0,0,0]),
        expression("irt_score",[[2,0]],[1,1,2,0,41,1,0,0,1,0,0]),
        expression("sem_ram",[[1,1,0]]),
        expression("sem_ram",[[0,-1,0]]),
        expression("sem_ram",[[0,1,1,0,0],[1,0,0,1,0]]),
        expression("sem_ram",[[0,0,-1,0,0],[0,0,0,1,0]]),
        expression("cfa_scores",[[1,2]],[0]),
        expression("degradation_path",[[0,1,0],[0,2,0]]),
        expression("degradation_predict",[[1,7,0]],[]),
        expression("destructive_predict",[[1,0,1]],[]),
        expression("cumulative_damage_predict",[[1,0,.5]],[]),
        "chem_irt_predict([[0.0]],NULL::DOUBLE[])",
        "chem_kinship_model(NULL::DOUBLE[][],[0.0])",
        "chem_irt_score([[0.0,NULL]],[1.0,1.0,2.0,0.0,41.0,1.0,0.0,0.0,1.0,0.0,0.0])",
        "chem_irt_predict([[0.0],[1.0,2.0]],[1.0,1.0,1.0,0.0,41.0,1.0,0.0,0.0])",
    ]
    bad=runner.query(" UNION ALL ".join(f"SELECT {i}AS id,({sql})IS NULL AS valid"for i,sql in enumerate(invalid))+" ORDER BY id;")
    assert len(bad)==len(invalid)and all(r["valid"]for r in bad),bad
    selected=runner.query("SELECT i,chem_irt_predict([[(i%3)::DOUBLE-1.0]],[1.0,3.0,1.0,0.0,41.0,2.0,0.0,0.2]) AS result FROM range(5000)t(i)WHERE i%997=0 ORDER BY i;")
    for row in selected:
        theta=row["i"]%3-1
        equal(strict_json(row["result"])["probabilities"],[[.2+.8/(1+math.exp(-2*theta))]],"selected IRT curve")
        runner.decoded+=1
    interleaved=runner.query("SELECT i,chem_irt_predict(CASE WHEN i%3=0 THEN NULL ELSE [[0.0]] END,[1.0,1.0,1.0,0.0,41.0,1.0,0.0,0.0]) AS result FROM range(8)t(i)WHERE i<>4 ORDER BY i;")
    for row in interleaved:
        if row["i"]%3==0:
            assert row["result"]is None
        else:
            equal(strict_json(row["result"])["probabilities"],[[.5]],"interleaved IRT")
            runner.decoded+=1
    print(json.dumps({"status":"PASS","layer":"native DuckDB SQL -> C++ -> Rust","specialized_catalog_examples":len(examples),"operations_with_numerical_fixtures":len(operations),"numeric_fields_checked":numeric_fields,"sql_saved_model_roundtrips":len(reuse),"static_ram_with_dynamic_instability_cases":1,"invalid_inputs_returning_sql_null":len(invalid),"selected_rows":len(selected),"interleaved_null_rows":len(interleaved),"native_sql_queries":runner.queries,"strict_json_results":runner.decoded},indent=2))


if __name__=="__main__":
    main()
