#!/usr/bin/env python3
"""Native DuckDB checks for time series, functional/text preprocessing and spatial models."""
import argparse,json
from pathlib import Path
from check_jmp_continuation import Runner,equal,expression,field,sql_list,strict_json
ROOT=Path(__file__).resolve().parents[1]
DEFAULT=['jmp-state-space-api.json','jmp-functional-preprocess-api.json','jmp-text-analysis-api.json','jmp-spatial-api.json']
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--duckdb',required=True);p.add_argument('--extension',required=True,type=Path);p.add_argument('--api',nargs='+',default=DEFAULT);p.add_argument('--expected',type=int,default=29);a=p.parse_args()
 runner=Runner(a.duckdb,a.extension);apis=[json.loads((ROOT/'tasks'/n).read_text())for n in a.api];ops=[o for g in apis for o in g['operations']];by={o['operation']:o for o in ops}
 assert len(ops)==len(by)==a.expected
 names=','.join("'"+o['sql_name']+"'"for o in ops)
 catalog=runner.query('SELECT function_name,description,examples,parameter_types FROM duckdb_functions() WHERE function_name IN('+names+');')
 assert len(catalog)==len(ops),(len(catalog),len(ops),'public SQL registration missing')
 assert all(o['description'] and o['examples'] and o['parameter_types']==['DOUBLE[][]','DOUBLE[]']for o in catalog)
 fields=0
 for op,result in zip(ops,runner.many([o['sql_example']for o in ops])):
  assert op['known_value'],op['operation']
  for key,value in op['known_value'].items():equal(field(result,key),value,op['operation']+'.'+key,5e-5);fields+=1
 reuse=[]
 def saved(fit,predict,query,expected,key='predict_parameters',prefix=()):
  if fit not in by:return
  params="json_extract(value,'$."+key+"')::DOUBLE[]"
  if prefix:params='list_concat('+sql_list(prefix)+','+params+')'
  sql='(WITH fit AS(SELECT '+by[fit]['sql_example']+' AS value)SELECT chem_'+predict+'('+sql_list(query)+','+params+')FROM fit)'
  reuse.append((fit,predict,sql,expected))
 saved('kalman_smoother','state_space_predict',[[0]],{'forecast_means':[1.4242424242424243]*2,'forecast_observation_variances':[.9390572390572391,1.1390572390572393]},'prediction_model',[2,.95])
 saved('text_weight','text_weight_predict',[[0,4,2]],{'matrix':[[.9542425094393249,0,.704365036222725]],'dictionary_refitted':False},'weight_parameters')
 saved('text_lsa','text_projection',[[1.5,1.5]],{'scores':[[.5,.5]],'preprocessing_refitted':False},'projection_parameters')
 saved('text_topic_svd','text_projection',[[1.5,1.5]],{'scores':[[.5,.75]],'preprocessing_refitted':False},'projection_parameters')
 saved('text_latent_class','text_latent_class_predict',[[1,0,1]],{'probabilities':[[1]]})
 saved('spatial_gls','spatial_predict',[[0],[10]],{'predictions':[5/3,3],'latent_variances':[7/9,3],'observation_variances':[16/9,4]})
 for fit,predict,sql,expected in reuse:
  result=runner.one(sql)
  for key,value in expected.items():equal(field(result,key),value,'saved '+fit+' -> '+predict+'.'+key,5e-5);fields+=1
 invalid=[
  expression('spatial_covariance',[[0],[1]],[3,1,0,1]),expression('spatial_gls',[[0,1],[0,2],[1,3]],[1,0,1,1,0,1]),expression('spatial_predict',[[0]],[]),
  expression('kalman_smoother',[[1,1]],[1,1,.95,0,1,1,-.1,1,.5]),expression('state_space_predict',[[0]],[]),
  expression('arima_forecast_intervals',[[1]],[.95,-1,0,0,0,1,1]),expression('cross_correlation',[[1,1],[1,1]],[1]),
  expression('text_weight',[[1,-1]],[]),expression('text_projection',[[1]],[]),expression('text_latent_class_predict',[[1,0]],[]),
 ]
 # These generic wrapper cases are different from each algorithm's domain rejection.
 invalid.extend([o['sql_name']+'(NULL::DOUBLE[][],'+sql_list(o['fixture']['args'])+')'for o in ops])
 valid=runner.query(' UNION ALL '.join('SELECT '+str(i)+' AS id,('+s+')IS NULL AS valid'for i,s in enumerate(invalid))+' ORDER BY id;')
 assert len(valid)==len(invalid) and all(v['valid']for v in valid),valid
 # Selected rows exercise nested-list dictionary vectors, interleaved validity and a constant parameter vector.
 selected=runner.query("SELECT i,chem_spatial_covariance([[0.0],[((i%3)+1)::DOUBLE]],[1.0,2.0,1.0,2.0]) AS result FROM range(5000)t(i)WHERE i%997=0 ORDER BY i;")
 import math
 for row in selected:equal(strict_json(row['result'])['covariance'][0][1],2*math.exp(-((row['i']%3)+1)/2),'selected covariance');runner.decoded+=1
 interleaved=runner.query("SELECT i,chem_text_projection(CASE WHEN i%3=0 THEN NULL ELSE [[1.5,1.5]] END,[1.0,2.0,2.0,3.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0])AS result FROM range(8)t(i)WHERE i<>4 ORDER BY i;")
 for row in interleaved:
  if row['i']%3==0:assert row['result']is None
  else:equal(strict_json(row['result'])['scores'],[[.5,.5]],'interleaved text projection');runner.decoded+=1
 print(json.dumps({'status':'PASS','layer':'native DuckDB SQL -> C++ -> Rust','catalog_examples':len(ops),'operations_with_numerical_fixtures':len(ops),'numeric_fields_checked':fields,'sql_saved_model_roundtrips':len(reuse),'invalid_inputs_returning_sql_null':len(invalid),'selected_rows':len(selected),'interleaved_null_rows':len(interleaved),'native_sql_queries':runner.queries,'strict_json_results':runner.decoded},indent=2))
if __name__=='__main__':main()
