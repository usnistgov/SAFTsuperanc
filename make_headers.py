import json
import textwrap
import ChebTools
import numpy as np

template_crit_head = """#pragma once

#include "ChebTools/ChebTools.h"
 
namespace PCSAFTSuperAncillary{ 
  Eigen::ArrayXd toarr(const std::vector<double>& e){ return Eigen::Map<const Eigen::ArrayXd>(&(e[0]), e.size()); }

  ChebTools::ChebyshevCollection cc_Ttilde{{
"""
template_crit_foot = """  }}; 
} /* namespace PCSAFTSuperAncillary */"""

def make_collection_contents(exps):
    s = ''
    for exp in exps:
        coef = str(exp['coef']).replace('[','{').replace(']','}')
        xmin, xmax = exp['xmin'], exp['xmax']
        line = f'    {{ toarr({coef:s}), {xmin}, {xmax} }},'
        s += line + '\n'
    s = s.strip(',')
    return s

def make_collection_contents_LV(exps):
    sL, sV = '',''
    for exp in exps:
        xmin, xmax = exp['xmin'], exp['xmax']

        coefL = str(exp['coefL']).replace('[','{').replace(']','}')
        line = f'  {{ std::vector<double>({coefL:s}), {xmin}, {xmax} }},'
        sL += line + '\n'

        coefV = str(exp['coefV']).replace('[','{').replace(']','}')
        line = f'  {{ std::vector<double>({coefV:s}), {xmin}, {xmax} }},'
        sV += line + '\n'

    sL = sL.strip(',')
    sV = sV.strip(',')
    return sL, sV

def write_crit():
    j = json.load(open('output/PCSAFT_crit_pts_expansions.json', 'r'))

    s = template_crit_head
    s += make_collection_contents(j['Ttilde'])
    s += template_crit_foot
    # print(s)
    with open('PCSAFTsuperancillary_crit.hpp', 'w') as fp:
        fp.write(s)

template_domain_head = """namespace PCSAFTSuperAncillary{ 
  auto toarr2(const std::initializer_list<double>& e){ return e; }
"""
template_domain_foot = """  CompleteWInterval domain{
    Wedges,
    {interval0,interval1,interval2,interval3,interval4,interval5,interval6,interval7,interval8,interval9}
  }; 
} /* namespace PCSAFTSuperAncillary */"""

def build_WInterval(N, wmin, wmax):
    s = rf'{{ // domain in w: [{wmin},{wmax}]' + '\n'
    wnodes = ChebTools.ChebyshevExpansion(np.linspace(0,1,17),wmin,wmax).get_nodes_realworld().tolist()
    coef = str(wnodes).replace('[','{').replace(']','}')
    s += '  // Wnodes\n'
    s += '  ' + coef + ',\n'
    s += '  { // Intervals (liq)\n'
    for w in wnodes:
        m = 1/w
        expansions = json.load(open(f'output/PCSAFT_VLE_m{m:0.12e}_expansions.json'))
        sL, sV = make_collection_contents_LV(expansions)
        s += textwrap.indent('{{\n' + sL + '\n}},\n', ' '*4) 
        # break
    s += '  },\n'
    s += '  { // Intervals (vapor)\n'
    for w in wnodes:
        m = 1/w
        expansions = json.load(open(f'output/PCSAFT_VLE_m{m:0.12e}_expansions.json'))
        sL, sV = make_collection_contents_LV(expansions) 
        s += textwrap.indent('{{\n' + sV + '\n}},\n', ' '*4) 
        # break
    s += '  },\n'
    s += '};\n'
    return s

def write_domain():
    # Get last Wedges file
    Wedges = json.load(open('output/Wedges_pass8.json'))['Wedges']

    s = template_domain_head
    s += '  // Wedges\n'
    s += '  const std::vector<double> Wedges = ' + str(Wedges).replace('[','{').replace(']','}') + ';\n\n'
    s += '  // Intervals\n'
    
    for i in range(len(Wedges)-1): # iterate over intervals in w
        interval = textwrap.indent(build_WInterval(16, Wedges[i], Wedges[i+1]), ' '*2)
        s += f'  const WInterval interval{i}     ' + interval + '\n'
        
    s += template_domain_foot

    # print(s)
    with open('domain.cpp', 'w') as fp:
        fp.write(s)

if __name__ == '__main__':
    write_crit()
    write_domain()