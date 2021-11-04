eqbm_cond_list = {
      'MUc' % marginal utility of consumption proxy
      'HHFOCc' % household FOC w.r.t consumption
      'HHFOCi' % household FOC w.r.t investment
      'HHFOCu' % household FOC w.r.t capital utilization
      'HHFOCk' % household FOC w.r.t capital
      'Euler' % intertemporal Euler equation
      'BondReturn' % longterm govt bond return
      'NAC' % No arbitrage condition
      'WageMarkup' % wage markup equation
      'Phillipsw' % wage Phillips curve
      'LOMk' % Law of motion for capital
      'Kutil' % capital utilization equation
      'MRTS' % marginal rate of technology substitution
      'PriceMarkup' % price markup equation
      'Phillipsp' % price Phillips curve
      'AggMktClr' % aggregate market clearing condition
      'AggProdFcn' % aggregate production function
      'MP' % monetary policy rule
      'TAXREV' % total tax revenue
      'PSDEF' % primary surplus definition
      'GBC' % government budget constraint
      'FPg' % fiscal policy rule: goverment spending
      'FPz' % fiscal policy rule: lump-sum transfer
      'MKVB' % market value of debt equation
      'ARgamma' % exogenous process for permenant growth
      'ARU' % exogenous process for int-pref
      'ARi' % exogenous process for marginal efficiency of investment
      'ARw' % exogenous process for wage markup
      'ARp' % exogenous process for price markup
      'ARg' % exogenous process for government spending
      'ARz' % exogenous process for lump-sum tax
      'ARM' % exogenous process for MP shock
      
      'EXgamma' % growth expectation
      'EXUc' % marginal utility of consumption expectation
      'EXi' % investment expectation
      'EXrk' % rental rate expectation
      'EXq' % Tobin's q expectation
      'EXlambda' % lagrange multiplier expectation
      'EXRL' % long term bond return expectation
      'EXinfl' % inflation expectation
      'EXw' % wage expectation
      
      'Lagy'
      'Lagivst'
      'Lagw'
      'Lagb'            
    };
n_eqbm_cond = length(eqbm_cond_list);
for eqbm_index = 1:n_eqbm_cond
    eval([char(eqbm_cond_list(eqbm_index)) ' = eqbm_index;']);
end