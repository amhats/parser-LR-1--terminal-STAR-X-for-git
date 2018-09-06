
#include "SMZP_Support.h"


/*****************
 * Alternating Array helpers
 *****************/
static void printAAZ(AltArrZ_t* aa) {
	for (int i = 0; i < AA_SIZE(aa); ++i) {
		gmp_fprintf(stderr, "%Zd*%llx + ", aa->elems[i].coef, aa->elems[i].degs);
	}
	fprintf(stderr, "\n");
}


int isExactlyEqual_AAZ(AltArrZ_t* a, AltArrZ_t* b) {
	if (a == NULL) {
		return (b == NULL);
	}
	if (b == NULL) {
		return 0;
	}

	if (a->size != b->size) {
		return 0;
	}

	for (int i = 0; i < a->size; ++i) {
		if (a->elems[i].degs != b->elems[i].degs) {
			return 0;
		}
		if (mpz_cmp(a->elems[i].coef, b->elems[i].coef) != 0) {
			return 0;
		}
	}

	return 1;
}

void expandNumVars_AAZ(AltArrZ_t* aa, int newNvar) {
	if (newNvar <= aa->nvar) {
		return;
	}
	if (aa->nvar == 0) {
		aa->nvar = newNvar;
		return;
	}

	degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);

	int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	int* __restrict__ newSizes = getExpOffsetArray(newNvar);

	degrees_t* maxExps = getMaxExpArray(newNvar);

	degrees_t degs;
	degrees_t curDeg;
	int diff = newNvar - aa->nvar;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < aa->nvar; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]); 
			if (curDeg > maxExps[j]) {
				fprintf(stderr, "SMZP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j]);  
		}
	}

	free(oldMasks);
	free(oldSizes);
	free(newSizes);

	aa->nvar = newNvar;
}

void expandNumVarsLeft_AAZ(AltArrZ_t* aa, int newNvar) {
	if (newNvar <= aa->nvar) {
		return;
	}
	if (aa->nvar == 0) {
		aa->nvar = newNvar;
		return;
	}

	degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);

	int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	int* __restrict__ newSizes = getExpOffsetArray(newNvar);

	degrees_t* maxExps = getMaxExpArray(newNvar);

	degrees_t degs;
	degrees_t curDeg;
	int diff = newNvar - aa->nvar;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < aa->nvar; ++j) {
			curDeg = GET_NTH_EXP(degs, oldMasks[j], oldSizes[j]); 
			if (curDeg > maxExps[j+diff]) {
				fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for expand at index %d; %lld > %lld.\n", j+diff, curDeg, maxExps[j+diff]);
				exit(1);
			}
			aa->elems[i].degs |= (curDeg << newSizes[j+diff]);  
		}
	}

	free(oldMasks);
	free(oldSizes);
	free(newSizes);

	aa->nvar = newNvar;
}

void shrinkNumVarsAtIdx_AAZ(AltArrZ_t* aa, int idx) {
	if (aa == NULL || aa->size < 1 || aa->nvar < 1) {
		return;
	}
	if (aa->nvar == 1) {
		for (int i = 0; i < aa->size; ++i) {
			aa->elems[i].degs = 0;
		}
		return;
	}

	int nvar = aa->nvar;
	int newNvar = aa->nvar - 1;

	AAZElem_t* elems = aa->elems;
	degrees_t* __restrict__ oldMasks = getExpMaskArray(aa->nvar);

	int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	int* __restrict__ newSizes = getExpOffsetArray(newNvar);

	degrees_t newDegs;
	degrees_t curDegs;
	degrees_t deg;
	int j;
	for (int i = 0; i < aa->size; ++i ) {
		curDegs = elems[i].degs;
		newDegs = 0;
		//iterate this current nvar, skipping j = idx, repacking exponents.
		//when j > idx, exponents are shifted to be at an index one less than 
		//originally.
		for (j = 0; j < idx; ++j) {
			deg = GET_NTH_EXP(curDegs, oldMasks[j], oldSizes[j]);
			newDegs |= (deg << newSizes[j]);
		}
		for (j = idx+1; j < nvar; ++j) {
			deg = GET_NTH_EXP(curDegs, oldMasks[j], oldSizes[j]);
			newDegs |= (deg << newSizes[j-1]);
		}

		elems[i].degs = newDegs;
	}

	free(oldMasks);
	free(oldSizes);
	free(newSizes);

	aa->nvar = newNvar;
}

void shrinkAndReorderVars_AAZ(AltArrZ_t* aa, int* varMap, int varmapSize) {
	if (aa == NULL || aa->size < 1 || aa->nvar < 1) {
		return;
	}

	if (varmapSize > aa->nvar) {
		return;
	} 

	int newNvar = 0;
	int needSort = 0;
	int maxSoFar = -1;
	for (int i = 0; i < varmapSize; ++i) {
		if (varMap[i] >= 0) {
			++newNvar;

			//the below checks for reordering, setting needSort if reorder occurs.
			if (varMap[i] < maxSoFar) {
				needSort = 1;
			} else {
				maxSoFar = varMap[i];
			}
		}
	}

	unsigned long long int* oldMasks = getExpMaskArray(aa->nvar);
	int* __restrict__ oldSizes = getExpOffsetArray(aa->nvar);
	int* __restrict__ newSizes = getExpOffsetArray(newNvar);
	degrees_t* __restrict__ maxExps = getMaxExpArray(newNvar);

	AAZElem_t* elems = aa->elems;
	degrees_t newDegs, oldDegs, curDeg;
	int j;
	for (int i = 0; i < aa->size; ++i) {
		oldDegs = elems[i].degs;
		newDegs = 0;
		for (j = 0; j < varmapSize; ++j) {
			if (varMap[j] < 0) {
				continue;
			}
			curDeg = GET_NTH_EXP(oldDegs, oldMasks[j], oldSizes[j]);
			if (curDeg > maxExps[varMap[j]]) {
				fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for shrink and reorder vars. At new index %d.", varMap[j]);
				exit(1);
			}
			newDegs |= (curDeg << newSizes[varMap[j]]);
		} 
		elems[i].degs = newDegs;
	}

	free(oldMasks);
	free(oldSizes);
	free(newSizes);
	free(maxExps);

	aa->nvar = newNvar;

	if (needSort) {
		mergeSortPolynomial_AAZ(aa);
	}
}

//     degs[varmap[i]] = curDegs[i];
void reorderVars_AAZ(AltArrZ_t* aa, int* varMap, int varmapSize) {
	int nvar = aa->nvar;

	degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	int* __restrict__ sizes = getExpOffsetArray(nvar);
	degrees_t* maxExps = getMaxExpArray(nvar);

	degrees_t degs;
	degrees_t curDeg;
	for (int i = 0; i < aa->size; ++i) {
		degs = aa->elems[i].degs;
		aa->elems[i].degs = 0;
		for (int j = 0; j < varmapSize; ++j) {
			curDeg = GET_NTH_EXP(degs, masks[j], sizes[j]);
			if (curDeg > maxExps[varMap[j]]) {
				fprintf(stderr, "SMQP ERROR: Overflow in exponent packing for reorder at index %d; %lld > %lld.\n", varMap[j], curDeg, maxExps[varMap[j]]);
				exit(1);	
			}
			aa->elems[i].degs |= (curDeg << sizes[varMap[j]]);
		}
	}

	mergeSortPolynomial_AAZ(aa);
}

degrees_t calculateMaxDegs_AAZ(AltArrZ_t* aa) {
	if (aa->nvar == 0) {
		return 0;
	}
	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	register int nvar = aa->nvar;
	degrees_t* __restrict__ masks = getExpMaskArray(nvar);
	degrees_t* __restrict__ maxList = (degrees_t*) calloc(nvar, sizeof(degrees_t));
	degrees_t max = 0ll;
	for (register int i = 0; i < size; ++i) {
		for (register int j = 0; j < nvar; ++j) {
			maxList[j] = (elems[i].degs & masks[j]) > (maxList[j]) ? (elems[i].degs & masks[j]) : maxList[j];  
		}
	}

	for (register int j = 0; j < nvar; ++j) {
		max |= maxList[j];
	}

	free(maxList);
	free(masks);

	return max;
}

AltArrZ_t* deepCopyPolynomial_AAZFromNode(Node* a, int nvar) {
	if (a == NULL) {
		return NULL;
	}
	polysize_t asize = numberOfTermsNode(a);
	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	AA_SIZE(newAA) = asize;
	newAA->alloc = asize;
	newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*newAA->alloc);
	
	int size = AA_SIZE(newAA);
	AAZElem_t* elems = newAA->elems;
	for (int i = 0; i < size && a != NULL; ++i) {
		mpz_init(elems[i].coef);
		mpz_set(elems[i].coef, mpq_numref(a->coef));
		elems[i].degs = a->degs;
		a = a->next;
	}	

	newAA->nvar = nvar;
	return newAA;
}

Node* deepCopyPolynomial_NodeFromAAZ(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}
	int asize = AA_SIZE(aa);
	
	AAZElem_t* elems = aa->elems;
	mpq_t qCoef;
	mpq_init(qCoef);
	mpz_set(mpq_numref(qCoef), elems[0].coef);
	Node* head = addTerm(NULL, elems[0].degs, qCoef);
	Node* tail = head;
	for (int i = 1; i < asize; ++i) {
		mpz_set(mpq_numref(qCoef), elems[i].coef);
		tail = addTerm(tail, elems[i].degs, qCoef);
	}	
	mpq_clear(qCoef);
	return head;
}

AltArr_t* deepCopyPolynomial_AAFromAAZ(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	AltArr_t* newAA = (AltArr_t*) malloc(sizeof(AltArr_t));
	AA_SIZE(newAA) = AA_SIZE(aa);
	newAA->alloc = aa->alloc;
	if (newAA->alloc > 0) {
		newAA->elems = (AAElem_t*) malloc(sizeof(AAElem_t)*newAA->alloc);
		int size = AA_SIZE(newAA);
		AAElem_t* elems = newAA->elems;
		AAZElem_t* oldelems = aa->elems;
		for (int i = 0; i < size; ++i) {
			mpq_init(elems[i].coef);
			mpz_set(mpq_numref(elems[i].coef), oldelems[i].coef);
			elems[i].degs = oldelems[i].degs;
		}
	} else {
		newAA->elems = NULL;
	}
	newAA->nvar = aa->nvar;

	return newAA;
}

AltArrZ_t* deepCopyPolynomial_AAZFromAA(AltArr_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	AA_SIZE(newAA) = AA_SIZE(aa);
	newAA->alloc = aa->alloc;
	if (newAA->alloc > 0) {
		newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*newAA->alloc);
		int size = AA_SIZE(newAA);
		AAZElem_t* elems = newAA->elems;
		AAElem_t* oldelems = aa->elems;
		for (int i = 0; i < size; ++i) {
			if (mpz_cmp_si(mpq_denref(oldelems[i].coef), 1l) != 0) {
				fprintf(stderr, "SMZP ERROR: Failed to convert a rational number polynomial to an integer one.\n");
				exit(1);
			}
			mpz_init(elems[i].coef);
			mpz_set(elems[i].coef, mpq_numref(oldelems[i].coef));
			elems[i].degs = oldelems[i].degs;
		}
	} else {
		newAA->elems = NULL;
	}
	newAA->nvar = aa->nvar;

	return newAA;

}

AltArrZDegList_t* deepCopyPolynomial_AAZDegListFromAA(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	int nvar = aa->nvar;
	unsigned long long int* masks = getExpMaskArray(nvar);
	int* sizes = getExpOffsetArray(nvar);

	AltArrZDegList_t* ret = makePolynomial_AAZDL(aa->size, nvar);
	AAZElem_DegList_t* rElems = ret->elems;
	AAZElem_t* elems = aa->elems;

	for (int i = 0; i < aa->size; ++i) {
		mpz_init(rElems[i].coef);
		mpz_set(rElems[i].coef, elems[i].coef);
		rElems[i].degs = malloc(sizeof(degree_t)*nvar);
		for (int j = 0; j < nvar; ++j) {
			rElems[i].degs[j] = GET_NTH_EXP(elems[i].degs, masks[j], sizes[j]);
		}
	}

	ret->size = aa->size;

	free(masks);
	free(sizes);

	return ret;
}

AltArrZ_t* deepCopyPolynomial_AAZ(AltArrZ_t* aa) {
	if (aa == NULL || AA_SIZE(aa) == 0) {
		return NULL;
	}

	AltArrZ_t* newAA = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
	AA_SIZE(newAA) = AA_SIZE(aa);
	newAA->alloc = aa->alloc;
	if (newAA->alloc > 0) {
		newAA->elems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*newAA->alloc);
		int size = AA_SIZE(newAA);
		AAZElem_t* elems = newAA->elems;
		AAZElem_t* oldelems = aa->elems;
		for (int i = 0; i < size; ++i) {
			mpz_init(elems[i].coef);
			mpz_set(elems[i].coef, oldelems[i].coef);
			elems[i].degs = oldelems[i].degs;
		}
	} else {
		newAA->elems = NULL;
	}
	newAA->nvar = aa->nvar;

	return newAA;
}

AltArrZ_t* sortPolynomial_AAZ(AltArrZ_t* aa) {
	//TODO not insertion sort.
	AAZElem_t* elems = aa->elems;
	int size = AA_SIZE(aa);

	degrees_t swapDegs;
	for (int i = 1; i < size; ++i) {
		for (int j = i; j > 0 && compareExponentVectors(elems[j-1].degs, elems[j].degs) < 0; --j) {
			mpz_swap(elems[j-1].coef, elems[j].coef);
			swapDegs = elems[j-1].degs;
			elems[j-1].degs = elems[j].degs;
			elems[j].degs = swapDegs;
		}
	}

	condensePolyomial_AAZ(aa);

	return aa;
}

static void mergeAAZElems(AAZElem_t* __restrict__ a, AAZElem_t* __restrict__ endA, AAZElem_t* __restrict__ b , AAZElem_t* __restrict__ endB, AAZElem_t* __restrict__ sorted) {
	int i = 0;
	while (a < endA && b < endB) {
		if (isGreaterExponentVectors(a->degs, b->degs)) {
			sorted[i] = *a; 
			++a;
		} else {
			sorted[i] = *b;
			++b;
		}
		++i;
	}

	while (a < endA) {
		sorted[i] = *a;
		++a;
		++i;
	}

	while (b < endB) {
		sorted[i] = *b;
		++b;
		++i;
	}
}

void mergeSortPolynomial_AAZ(AltArrZ_t* aa) {
	if (aa->size < 9) {
		sortPolynomial_AAZ(aa);
	}

	int size = aa->size;
	AAZElem_t* tempElems = (AAZElem_t*) malloc(sizeof(AAZElem_t)*size);
	AAZElem_t* elems = aa->elems;
	int end1 = 0;
	int end2 = 0;
	for (int window = 1; window < size; window <<= 1) {
		//merge elems[i]:elems[i+window] with elems[i+window]:elems[i+2*window]
		for (int i = 0; i < size; i += 2*window) {
			end1 = i + window < size ? i + window : size;
			end2 = i + 2*window < size ? i + 2*window : size;
			mergeAAZElems(elems+i, elems+end1, elems+end1, elems+end2, tempElems+i);
		}

		AAZElem_t* temp = tempElems;
		tempElems = elems;
		elems = temp;
	}

	aa->elems = elems;
}

void condensePolyomial_AAZ(AltArrZ_t* aa) {
	if (AA_SIZE(aa) < 1) {
		return;
	}

	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	int insertIdx = 0;
	int compareIdx = 1;
	while (compareIdx < size) {
		if (compareExponentVectors(elems[insertIdx].degs, elems[compareIdx].degs) == 0) {
			mpz_add(elems[insertIdx].coef, elems[insertIdx].coef, elems[compareIdx].coef);
		} else if(compareIdx - insertIdx > 1) {
			++insertIdx;
			elems[insertIdx].degs = elems[compareIdx].degs;
			mpz_swap(elems[insertIdx].coef, elems[compareIdx].coef);
		} else {
			++insertIdx;
		}
		++compareIdx;		
	}

	++insertIdx;
	for (int i = insertIdx; i < size; ++i) {
		mpz_clear(elems[i].coef);
	}
	AA_SIZE(aa) = insertIdx;
}

void negatePolynomial_AAZ(AltArrZ_t* aa) {
	int size = AA_SIZE(aa);
	AAZElem_t* elems = aa->elems;
	for (int i = 0; i < size; ++i) {
		mpz_neg(elems[i].coef, elems[i].coef);
	}
}

void evalPolyToVal_AAZ(AltArrZ_t* aa, mpz_t* vals, int nvar, mpz_t res) {
	if (aa == NULL || nvar != aa->nvar) {
		mpz_set_ui(res, 0ul);
		return;
	}

	if (nvar == 0 || aa->nvar == 0) {
		mpz_set(res, aa->elems->coef);
		return;
	}

	int* sizes = getExpOffsetArray(nvar);
	unsigned long long int* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AAZ(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	mpz_set_ui(res, 0ul);
	mpz_t acc; 
	mpz_init(acc);

	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpz_set(acc, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			mpz_mul(acc, acc, valList[j][deg]);
		}
		mpz_add(res, res, acc);
	}

	mpz_clear(acc);

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	free(masks);
	free(sizes);
}

AltArrZ_t* evaluatePoly_AAZ(AltArrZ_t* aa, int* active, mpz_t* vals, int nvar) {
	if (aa == NULL) {
		return NULL;
	}

	if (nvar == 0 || aa->nvar == 0) {
		return deepCopyPolynomial_AAZ(aa);
	}

	int i;
	int newNvar = nvar;
	for (i = 0; i < nvar; ++i) {
		if (active[i]) {
			--newNvar;
		}
	}

	if (newNvar == 0) {
		mpz_t res;
		mpz_init(res);
		evalPolyToVal_AAZ(aa, vals, nvar, res);
		AltArrZ_t* ret = makeConstPolynomial_AAZ(1, 0, res);
		mpz_clear(res);
		return ret;
	} 

	int* sizes = getExpOffsetArray(nvar);
	unsigned long long int* masks = getExpMaskArray(nvar);

	mpz_t* valList[nvar];
	int valListSize[nvar];

	degrees_t maxDegs = calculateMaxDegs_AAZ(aa);
	for (int j = 0; j < nvar; ++j) {
		degree_t deg = GET_NTH_EXP(maxDegs, masks[j], sizes[j]);
		if (!active[j]) {
			valList[j] = NULL;
			valListSize[j] = 0;
			continue;
		}
		valList[j] = (mpz_t*) malloc(sizeof(mpz_t)*(deg+1));
		valListSize[j] = deg+1;

		mpz_init(valList[j][0]);
		mpz_set_ui(valList[j][0], 1ul);
		for (int k = 1; k < deg+1; ++k) {
			mpz_init(valList[j][k]);
			mpz_mul(valList[j][k], valList[j][k-1], vals[j]);
		}
	}

	AltArrZ_t* res = makePolynomial_AAZ(aa->size, newNvar);
	AAZElem_t* resElems = res->elems;
	res->size = aa->size;

	int* newSizes = getExpOffsetArray(newNvar);
	int size = aa->size;
	int k = 0;
	degrees_t newDegs = 0;
	for (int i = 0; i < size; ++i) {
		mpz_init(resElems[i].coef);
		mpz_set(resElems[i].coef, aa->elems[i].coef);
		for (int j = 0; j < nvar; ++j) {
			degree_t deg = GET_NTH_EXP(aa->elems[i].degs, masks[j], sizes[j]);
			if (valList[j] == NULL) {
				newDegs |= deg << newSizes[k];
				++k;
			} else {
				mpz_mul(resElems[i].coef, resElems[i].coef, valList[j][deg]);
			}
		}
		resElems[i].degs = newDegs;
		k = 0;
		newDegs = 0;
	}

	for (int j = 0; j < nvar; ++j) {
		for (int k = 0; k < valListSize[j]; ++k) {
			mpz_clear(valList[j][k]);
		}
		free(valList[j]);
	}

	mergeSortPolynomial_AAZ(res);

	return res;
}


AltArrZ_t* convertFromAAZElemToAAZ (AAZElem_t* coef, int coefSize, int nvar)
{
    if (coef == NULL || coefSize == 0) {
      return NULL;
    }

    AltArrZ_t* poly = makePolynomial_AAZ (coefSize, nvar);
    poly->size = coefSize;
    poly->elems = coef;
    return poly;
}

AltArrZ_t* swappingExponents_AAZ (AltArrZ_t* aa, int idx1, int idx2)
{
  if (aa == NULL)
      return NULL;

  
  AltArrZ_t* cPoly = deepCopyPolynomial_AAZ (aa);
  if (idx1 == idx2){
      return cPoly;
  }

  int varMap[cPoly->nvar];
  for (int i = 0; i < cPoly->nvar; ++i){
    varMap[i] = i;
  }
  varMap[idx1] = idx2;
  varMap[idx2] = idx1;

  reorderVars_AAZ (cPoly, varMap, cPoly->nvar);
  return cPoly;
}


int mainLeadingDegree_AAZ (AltArrZ_t* aa)
{
    if (aa == NULL || aa->size == 0){
	return 0;
    }
    
    int mvarDegOffset = getMVarExpOffset (aa->nvar);
    unsigned long long int* masks = getExpMaskArray (aa->nvar);
    unsigned long long int mvarMask = masks[0];
    free(masks);
    
    return ((aa->elems[0].degs & mvarMask) >> mvarDegOffset);
}


AltArrZ_t* mainLShiftPolynomial_AAZ (AltArrZ_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
	return NULL;
    }
    if (n < 1){
	return deepCopyPolynomial_AAZ (aa);
    }
    
    int nvar = aa->nvar;
    int mvarDegOffset = getMVarExpOffset (nvar);
    
    AltArrZ_t* shifted = (AltArrZ_t*) malloc (sizeof (AltArrZ_t));
    shifted->nvar = aa->nvar;
    shifted->size = aa->size;
    shifted->alloc = aa->alloc;
    shifted->elems = (AAZElem_t*) malloc (sizeof (AAZElem_t)*(aa->size));
    AAZElem_t* selems = shifted->elems;
    AAZElem_t* aelems = aa->elems;

    for (int i = 0; i < shifted->size; ++i){
	mpz_init (selems[i].coef);
	mpz_set (selems[i].coef, aelems[i].coef);
	selems[i].degs = aelems[i].degs;
	selems[i].degs += ((degrees_t) n << mvarDegOffset);
    }   
    return shifted;
}


AltArrZ_t* mainLShiftPolynomial_AAZ_inp (AltArrZ_t* aa, int n)
{
    if (aa == NULL || aa->size == 0){
	return NULL;
    }
    if (n < 1){
	return aa;
    }
    
    int nvar = aa->nvar;
    int mvarDegOffset = getMVarExpOffset (nvar);
    
    AAZElem_t* aelems = aa->elems;

    for (int i = 0; i < aa->size; ++i){
	aelems[i].degs += ((degrees_t) n << mvarDegOffset);
    }   
    return aa;
}

AltArrZ_t* leadingTerm_AAZ (AltArrZ_t* aa, int nvar)
{
	if (aa == NULL || aa->size == 0){
		return NULL;
	}

	AltArrZ_t* lt = makeConstPolynomial_AAZ (1, aa->nvar, aa->elems->coef);
	lt->elems->degs = aa->elems->degs;
	
	return lt;
}

int leadingVariable_AAZ (AltArrZ_t* aa)
{
	if (aa == NULL || aa->size == 0){
		return -2;
	}

	unsigned long long int* masks = getExpMaskArray(aa->nvar);
	for (int i = 0; i < aa->nvar; ++i){
		if ((aa->elems[0].degs & masks[i]) != 0) {
			return i;
		}
	}

	return -1;
}


AltArrZ_t* mainLeadingCoefficient_AAZ (AltArrZ_t* aa)
{
    if (aa == NULL || aa->size == 0){
	return NULL;
    }
    
    int mvarDegOffset = getMVarExpOffset (aa->nvar);
    unsigned long long int* masks = getExpMaskArray (aa->nvar);
    unsigned long long int mvarMask = masks[0];
    free(masks);
    
    AAZElem_t* elems = aa->elems;
    degree_t mvarDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
    degree_t curDeg = mvarDeg;
    
    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);
    
    for (int i = 0; (mvarDeg == curDeg) && i < AA_SIZE (aa); ++i){
	
	if (polyElemsSize + 1 > polyElemsAlloc){
	    polyElemsAlloc += 10;
	    polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t) * polyElemsAlloc);
	}

	mpz_init (polyElems[i].coef);
	mpz_set (polyElems[i].coef, elems[i].coef);
	polyElems[i].degs = (elems[i].degs & (~mvarMask));
	++polyElemsSize;
		
	if (i+1 < AA_SIZE(aa)){
	    curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
	} else {
	    curDeg = -1;
	}
    }

    AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
    poly->alloc = polyElemsAlloc;
    poly->size = polyElemsSize;
    poly->elems = polyElems;
    poly->nvar = aa->nvar;

    return poly;
}

AltArrZ_t* mainCoefficientAtIdx_AAZ (AltArrZ_t* aa, int e)
{
    if (aa == NULL || aa->size == 0){
	return NULL;
    }
    if (e < 0){
	return NULL;
    }
    
    int mvarDegOffset = getMVarExpOffset (aa->nvar);
    unsigned long long int* masks = getExpMaskArray (aa->nvar);
    unsigned long long int mvarMask = masks[0];
    free(masks);

    AAZElem_t* elems = aa->elems;
    degree_t mvarDeg = e;
    degree_t curDeg = (elems[0].degs & mvarMask) >> mvarDegOffset;
    int begin = -1;
    
    if (mvarDeg > curDeg){
        return NULL;
    }
    
    int polyElemsAlloc = 10;
    int polyElemsSize = 0;
    AAZElem_t* polyElems = (AAZElem_t*) malloc (sizeof(AAZElem_t) * polyElemsAlloc);
    
    for (int i = 0; i < AA_SIZE (aa); ++i){
	if (mvarDeg == curDeg){
	    if (polyElemsSize + 1 > polyElemsAlloc) {
		polyElemsAlloc += 10;
		polyElems = (AAZElem_t*) realloc (polyElems, sizeof (AAZElem_t)*polyElemsAlloc);
	    }
	    
	    mpz_init (polyElems[polyElemsSize].coef);
	    mpz_set (polyElems[polyElemsSize].coef, elems[i].coef);
	    polyElems[polyElemsSize].degs = elems[i].degs & (~mvarMask);
	    ++polyElemsSize;
	    
	    if (begin == -1) {
		begin = i;
	    }
	}
	
	if (begin != -1 && mvarDeg != curDeg){
	    break;
	}
	
	if (i+1 < AA_SIZE(aa)){
	    curDeg = (elems[i+1].degs & mvarMask) >> mvarDegOffset;
	} else {
	    curDeg = -1;
	}
    }
    
    AltArrZ_t* poly = (AltArrZ_t*) malloc(sizeof(AltArrZ_t));
    poly->alloc = polyElemsAlloc;
    poly->size = polyElemsSize;
    poly->elems = polyElems;
    poly->nvar = aa->nvar;
    
    return poly;
}

AltArrZ_t* maxPolynomials_AAZ (AltArrZ_t* a, AltArrZ_t* b)
{
  if (a == NULL || a->size == 0){
      if (b == NULL || b->size == 0){
          return NULL;
      }
      return deepCopyPolynomial_AAZ (b);
  }
  if (b == NULL || b->size == 0){
      return deepCopyPolynomial_AAZ (a);
  }

  int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

  if (cmp > 0){
      return deepCopyPolynomial_AAZ (a);
  } else if (cmp < 0){
      return deepCopyPolynomial_AAZ (b);
  } else {
      if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
          return deepCopyPolynomial_AAZ (a);
      } else {
          return deepCopyPolynomial_AAZ (b);
      }
  }
  return NULL;
}

AltArrZ_t* maxPolynomials_AAZ_inp (AltArrZ_t* a, AltArrZ_t* b)
{
  if (a == NULL || a->size == 0){
      if (b == NULL || b->size == 0){
          return NULL;
      }
      return deepCopyPolynomial_AAZ (b);
  }
  if (b == NULL || b->size == 0){
      return a;
  }

  int cmp = compareExponentVectors(a->elems[0].degs, b->elems[0].degs);

  if (cmp > 0){
      return a;
  } else if (cmp < 0){
      return deepCopyPolynomial_AAZ (b);
  } else {
      if (mpz_cmp(a->elems[0].coef, b->elems[0].coef) > 0) {
          return a;
      } else {
          return deepCopyPolynomial_AAZ (b);
      }
  }
  return NULL;
}



/*****************
 * SMZP Addition
 *****************/

AltArrZ_t* addPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b); 
	
	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);


	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_add(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);

	return c;
}


AltArrZ_t* addPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b); 
	
	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
			mpz_add(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			cElems[k] = aElems[i];
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		k += (asize - i);
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ(c, k);

	return c;
}

AltArrZ_t* subPolynomials_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b); 
	
	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);

	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	register int k = 0;
	register int i = 0;
	register int j = 0;

	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			mpz_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, aElems[i].coef, bElems[j].coef);
			cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			mpz_init(cElems[k].coef);
			mpz_set(cElems[k].coef, aElems[i].coef);
			cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	while(i < asize) {
		mpz_init(cElems[k].coef);
		mpz_set(cElems[k].coef, aElems[i].coef);
		cElems[k].degs = aElems[i].degs;
		++k;
		++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);

	return c;
}

AltArrZ_t* subPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	if (a == NULL && b == NULL) {
		return NULL;
	}
	register int asize = a == NULL ? 0 : AA_SIZE(a);
	register int bsize = b == NULL ? 0 : AA_SIZE(b); 
	
	AltArrZ_t* c = makePolynomial_AAZ(asize + bsize, nvar);
	
	AAZElem_t* aElems = a == NULL ? NULL : a->elems;
	AAZElem_t* bElems = b == NULL ? NULL : b->elems;
	AAZElem_t* cElems = c->elems;

	// ratNum_t ccoef;
	// mpq_init(ccoef);

	// register degree_t cmp;
	register int k = 0;
	register int i = 0;
	register int j = 0;
	
	while(i < asize && j < bsize) {
		if (isLessExponentVectors(aElems[i].degs, bElems[j].degs)) {
			//a < b
			mpz_init(cElems[k].coef);
			mpz_neg(cElems[k].coef, bElems[j].coef);
			cElems[k].degs = bElems[j].degs;
			++k;
			++j;
		} else if (isEqualExponentVectors(aElems[i].degs, bElems[j].degs)) {
			// a==b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			mpz_sub(cElems[k].coef, cElems[k].coef, bElems[j].coef);
			// cElems[k].degs = aElems[i].degs;
			if (mpz_sgn(cElems[k].coef) == 0) {
				mpz_clear(cElems[k].coef);
				// c[k] is a[k], so the above clears both.
				// mpq_clear(aElems[i].coef);
				--k;
			}
			++k;
			++i;
			++j;
		} else {
			//a > b
			cElems[k] = aElems[i];
			// mpq_init(cElems[k].coef);
			// mpq_set(cElems[k].coef, aElems[i].coef);
			// cElems[k].degs = aElems[i].degs;
			++k;
			++i;
		}
	}

	if(i < asize) {
		memcpy(cElems + k, aElems + i, sizeof(AAZElem_t)*(asize - i));
		// mpq_init(cElems[k].coef);
		// mpq_set(cElems[k].coef, aElems[i].coef);
		// cElems[k].degs = aElems[i].degs;
		k += (asize - i);
		// ++k;
		// ++i;
	}
	while(j < bsize) {
		mpz_init(cElems[k].coef);
		mpz_neg(cElems[k].coef, bElems[j].coef);
		cElems[k].degs = bElems[j].degs;
		++k;
		++j;
	}

	//free the arrays but do NOT free the underlying gmp data, used by c now.
	free(aElems);
	free(a);

	AA_SIZE(c) = k;

	resizePolynomial_AAZ(c, k);

	return c;
}



/*****************
 * SMQP Multiplication & Helpers
 *****************/

ProductHeap_AAZ* prodheapInit_AAZ(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	h->elements = (ProductHeapElem_AAZ*) malloc(sizeof(ProductHeapElem_AAZ)*AA_SIZE(a));
	h->elements[0].chain = prodheapMakeChain_AAZ(0, 0, NULL);
	addExponentVectors(a->elems->degs, b->elems->degs, h->elements->degs);
	h->heapSize = 1;
	h->maxHeapSize = AA_SIZE(a);
	h->nvar = nvar;
	return h;
}

/**
 * Insert a new element, elem, into the product heap, h, chaining as possible.
 */
void prodheapInsert_AAZ(ProductHeap_AAZ* h, ProductHeapChain_AAZ* chain, register degrees_t degs) {
	register int s = h->heapSize;
	ProductHeapElem_AAZ* elems = h->elements;

	if (s == 0) {
		elems[0].degs = degs;
		elems[0].chain = chain;
		h->heapSize = 1;
		return;		
	}
	
	//first check if we can chain off the root
	if (isEqualExponentVectors(elems[0].degs, degs)) {
		chain->next = elems[0].chain;
		elems[0].chain = chain;
		return;
	} 

	//otherwise, we must search the heap to find the new product's insertion point
	//note that since we are looking for chains we cannot use the simple swim method
	//we sort of fake the swimming, looking for a chain to be made or the eventual
	//place where the swim would stop. At this point, we insert the new elem
	//in that spot, and "push" the entire path we took down a level. Assuming 
	//that we insert e and it ends up at the root, we push down the 'x' path
	//                                      //
	//     x     --->    e                  //
	//    / \           / \                 //     
	//   x   o         x   o
	//                /
 	//               x

	register int i = (s-1) >> 1; //i is parent 
	register int j = s;       //j is current insertion point
	register long long unsigned int path = 1;
	while (j > 0) {
		if (isEqualExponentVectors(elems[i].degs, degs)) {
			chain->next = elems[i].chain;
			elems[i].chain = chain;
			return;
		} else if (isLessExponentVectors(elems[i].degs, degs)) {
			path <<= 1;
			if (!(j & 1)) {
				//set the trailing bit to 1 to distinguish left/right of path
				path += 1; 
			}
			j = i;
			i = (i-1) >> 1;
		} else { //cmp > 0
			break;
		}
	}

	//then j is now the place we need to insert elem;
	//do so, and then push all others down the path, inserting the last
	//as the new element in elems[s];
	ProductHeapElem_AAZ temp;
	ProductHeapElem_AAZ elem = {degs, chain};

	//TODO use i index again to swap between elements rather than use a second temp elem.
	while (j <= s) {
		temp = elems[j];
		elems[j] = elem;
		elem = temp;
		j = (j << 1) + 1 + (path & 1);
		path >>= 1;
	}
	++(h->heapSize);
}

/**
 * Extract the maximal heap element (chain) from the heap and return it.
 * Automatic insertion of the next element, a_i * b_j+1, is not done.
 * This allows the multiplication to limit the number of entries in the heap.
 * returns NULL if no such element in the heap.
 */
ProductHeapChain_AAZ* prodheapRemoveMax_AAZ(ProductHeap_AAZ* h) {
	ProductHeapElem_AAZ* elems = h->elements;
	ProductHeapChain_AAZ* maxElem = elems[0].chain;
	register int i = 0;
	register int j = 1;
	register int s = --(h->heapSize);
	
	//promote largest children
	while (j < s) {
		if (j+1 < s && isLessExponentVectors(elems[j].degs, elems[j+1].degs)) {
			++j;
		}
		elems[i] = elems[j];
		i = j;
		j = (j << 1) + 1;
	}
	//now place last element into i and swim up to make tree complete 
	j = (i-1) >> 1;
	while(i > 0) {
		if (isLessExponentVectors(elems[s].degs, elems[j].degs)) {
			break;
		}
		elems[i] = elems[j];
		i = j;
		j = (j-1) >> 1;
	}
	elems[i] = elems[s]; 

	return maxElem;
}

/**
 * Multiply two polynomials given their head Nodes, a and b.
 * This algorithm makes use of heaps as an efficient search data structure. 
 * It is assumed that both a and b have compatible exponent vectors.
 * 
 * nvar: number of elements in the exponent vectors.
 *
 * returns a pointer to the head Node of the product polynomial.
 */
AltArrZ_t* multiplyPolynomials_AAZ(AltArrZ_t* __restrict__ a, AltArrZ_t* __restrict__ b, int nvar) {
	if (a == NULL || a->size == 0 || b == NULL || b->size == 0) {
		return NULL;
	}

	degrees_t aMax = calculateMaxDegs_AAZ(a);
	degrees_t bMax = calculateMaxDegs_AAZ(b);

	checkValidMonomialMult(aMax, bMax, nvar);

	// reorder to obtain smaller as a. 
	if (b->size < a->size) {
		AltArrZ_t* temp = a;
		a = b;
		b = temp;
	}

	ProductHeap_AAZ* h = prodheapInit_AAZ(a,b,nvar);
	// mpz_t ccoef;
	// mpz_init(ccoef);
	
	//TODO smarter allocation here? dynamic reallocating? 
	AltArrZ_t* c = makePolynomial_AAZ(AA_SIZE(a)*AA_SIZE(b), nvar);
	
	//k is c, i is a, j is b.
	register int k = 0;
	// register int i = 0;
	// register int j = 0;

	AAZElem_t* __restrict__ cElems = c->elems;
	AAZElem_t* __restrict__ bElems = b->elems;

	AAZElem_t* aElems = a->elems;
	register int lastA = AA_SIZE(a) - 1;
	register int lastB = AA_SIZE(b) - 1; 
	register int firstB = 0; 	

	ProductHeapChain_AAZ* maxElem = NULL;
	ProductHeapChain_AAZ* nextMaxElem = NULL;
	degrees_t* nextDegs;
	while ( (nextDegs = prodheapPeek_AAZ(h)) != NULL) {
		//cache since, on RemoveMax, pointer is invalidated.
		cElems[k].degs = *nextDegs;
		mpz_init(cElems[k].coef);

		while (nextDegs != NULL && isEqualExponentVectors(cElems[k].degs, *nextDegs)) {
			//we will extract and accumulate the coefficents 
			//oldMaxElem and maxElem are both chains. We must merge both chains.
			//we do this by taking the head of maxElem and, as necessary, push it
			//to the head of the oldMaxElem chain
			ProductHeapChain_AAZ* oldMaxElem = maxElem;
			maxElem = prodheapRemoveMax_AAZ(h);
			while (maxElem != NULL) {

				mpz_addmul(cElems[k].coef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_mul(ccoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
				// mpz_add(cElems[k].coef, cElems[k].coef, ccoef);

				//If we extracted a_i*b_1 we need to insert a_(i+1)*b_1;
				if (maxElem->b == firstB && maxElem->a_i != lastA) {
					oldMaxElem = prodheapMakeChain_AAZ((maxElem->a_i)+1, firstB, oldMaxElem);
				}

				//cache next before freeing or overwriting 
				nextMaxElem = maxElem->next;

				//If the extracted term has another product in the stream, 
				//update the product and push onto the oldMaxElem chain
				if(maxElem->b != lastB) {
					++(maxElem->b);
					maxElem->next = oldMaxElem;
					oldMaxElem = maxElem;
				} else {
					//we are done with the maxElem ProductHeapChain
					maxElem->next = NULL;
					prodheapFreeChain_AAZ(maxElem);
				}

				maxElem = nextMaxElem;
			}

			//reset head of maxElem list
			maxElem = oldMaxElem;	

			nextDegs = prodheapPeek_AAZ(h);		
		}

		//Commit new term to the product.
		if (mpz_sgn(cElems[k].coef) != 0) {
			++k;
		} else {
			//reset accumulator variables and do not increment k.
			mpz_clear(cElems[k].coef);
		}
		
		//Insert all successors of previously extracted products
		while(maxElem != NULL) {
			//clear maxElem->next before inserting
			nextMaxElem = maxElem->next;
			maxElem->next = NULL;
			prodheapInsert_AAZ(h, maxElem, aElems[maxElem->a_i].degs + bElems[maxElem->b].degs);
			maxElem = nextMaxElem;
		}
	}

	// mpz_clear(ccoef);
	prodheapFree_AAZ(h);

	AA_SIZE(c) = k;
	resizePolynomial_AAZ(c, k);

	return c;
}

AltArrZ_t* multiplyPolynomials_AAZ_inp(AltArrZ_t* a, AltArrZ_t* b, int nvar) {
	//TODO actually do this in-place.
	AltArrZ_t* prod = multiplyPolynomials_AAZ(a,b,nvar);
	freePolynomial_AAZ(a);
	return prod;
}



/*****************
 * Polynomial exponentiation 
 *****************/


AltArrZ_t* exponentiatePoly_AAZ(AltArrZ_t* a, unsigned int n, int nvar) {
	if (n == 0) {
		AltArrZ_t* ret = makePolynomial_AAZ(1, nvar);
		mpz_init(ret->elems->coef);
		mpz_set_ui(ret->elems->coef, 1ul);
		ret->elems->degs = 0;
		ret->size = 1;
		return ret;
	} else if (n == 1) {
		return deepCopyPolynomial_AAZ(a);
	}

	AltArrZ_t* r = NULL;
	AltArrZ_t* b = deepCopyPolynomial_AAZ(a);
	while (n > 1) {
		if (n & 1) {
			r = (r == NULL) ? deepCopyPolynomial_AAZ(b) : multiplyPolynomials_AAZ_inp(r, b, nvar); 
		}
		b = multiplyPolynomials_AAZ_inp(b, b, nvar);
		n >>= 1;
	}
	r = (r == NULL) ? deepCopyPolynomial_AAZ(b) : multiplyPolynomials_AAZ_inp(r, b, nvar);

	freePolynomial_AAZ(b);
	return r;
}



/*****************
 * Polynomial division
 *****************/

/**
 * Extract a product term from the heap.
 * This product term is a_i*b_j for some i and j.
 * If the term b_(j+1) exists then the heap is updated by inserting a_i*b_(j+1).
 * This process continues as long as the next element in the heap has the same
 * product degree.
 */
void divisionGetNextTerm_AAZ(ProductHeap_AAZ* h, const AAZElem_t* __restrict__ aElems, const AAZElem_t* __restrict__ bElems, mpz_t* retCoef) {
	if (h->heapSize == 0) {
		return;
	}

	int lastB = h->lastB;

	ProductHeapChain_AAZ* insertChain = NULL;
	ProductHeapChain_AAZ* maxElem, *nextMaxElem;

	// mpz_t prodCoef;
	// mpz_init(prodCoef);
	degrees_t* nextDegs = prodheapPeek_AAZ(h);
	register degrees_t maxDegs = *nextDegs;

	while ( nextDegs != NULL && isEqualExponentVectors(maxDegs, *nextDegs)) {
		maxElem = prodheapRemoveMax_AAZ(h);
				
		while (maxElem != NULL) {
			nextMaxElem = maxElem->next;

			mpz_addmul(*retCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			// mpz_mul(prodCoef, aElems[maxElem->a_i].coef, bElems[maxElem->b].coef);
			// mpz_add(*retCoef, *retCoef, prodCoef);
			if (maxElem->b != lastB) {
				++(maxElem->b);
				maxElem->next = insertChain;
				insertChain = maxElem;
			} else {
				maxElem->next = NULL;
				prodheapFreeChain_AAZ(maxElem);
			}

			maxElem = nextMaxElem;
		}

		nextDegs = prodheapPeek_AAZ(h);
	}

	while(insertChain != NULL) {
		maxElem = insertChain->next;
		insertChain->next = NULL;
		prodheapInsert_AAZ(h,insertChain, aElems[insertChain->a_i].degs + bElems[insertChain->b].degs);
		insertChain = maxElem;
	}

	// mpz_clear(prodCoef);
}

void divideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int nvar) {
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	if (nvar == 0) {
		if(mpz_divisible_p(c->elems->coef, b->elems->coef)) {
			*res_a = makeConstPolynomial_AAZ(1, 0, c->elems->coef);
			mpz_divexact((*res_a)->elems->coef, (*res_a)->elems->coef, b->elems->coef);
			*res_r = NULL; 
		} else {
			*res_a = NULL;
			*res_r = makeConstPolynomial_AAZ(1, 0, c->elems->coef);
		}
		return;
	}


	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0; 

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar) && mpz_divisible_p(k->coef, bElems->coef)) {
			mpz_divexact(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = k->degs;
			++j;
			++(curR);
			mpz_init(curR->coef);
		}
		++k;
	}

	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i; 
	AA_SIZE(r) = j; 
	a->alloc = maxSize;
	r->alloc = maxSize;

	*res_a = a;
	*res_r = r;
}

/** 
 * Given two polynomials, c and b, find their quotient and remainder such that
 * c = b*a + r. The quotient a is returned in res_a, and the remainder r in res_r 
 * Based on Stephen Johnson's "Sparse Polynomial Arithmetic".
 */
void dividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, register int nvar) {

	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AAZ(0, nvar);
		*res_r = NULL;
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		divideBySingleTerm_AAZ(c, b, res_a, res_r, nvar);
		return;
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0; 

	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	
	register degrees_t beta = b->elems->degs;


	//init a with lt(c)/lt(b);
	DivTest_ptr divTest = getMonomialDivideTestFuncPtr(nvar);
	while (k != lenK && !( (*divTest)(k->degs, beta) && mpz_divisible_p(k->coef, b->elems->coef)) ) {
		mpz_set(curR->coef, k->coef);
		curR->degs = k->degs;
		++j;
		if (j >= maxRSize) {
			maxRSize <<= 1;
			r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
			curR = r->elems + j - 1; 	
		}
		++(curR);
		mpz_init(curR->coef);
		++k;
	}

	if (k == lenK) {
		//no division to do at all!
		mpz_clear(curA->coef);
		mpz_clear(curR->coef);

		AA_SIZE(a) = i;
		AA_SIZE(r) = j;
		a->alloc = maxASize;
		r->alloc = maxRSize;

		*res_a = a;
		*res_r = r;
		return; 
	}

	subtractExponentVectors(k->degs, beta, curA->degs);
	//assuming b is monic!
	// mpz_set(curA->coef, k->coef);
	mpz_divexact(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h); 
	register degrees_t eps;
	register cmpExp_t cmp;
	while(delta != NULL || k != lenK) {
		
		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if ((*divTest)(eps, beta) && mpz_divisible_p(curA->coef, b->elems->coef)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a  	
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {

			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				curR = r->elems + j - 1; 	
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	*res_a = a;
	*res_r = r;
} 

void exactDividePolynomials_AAZ (AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, register int nvar)
{
    
    if (b == NULL || AA_SIZE(b) == 0) {
	//division by zero
	fprintf (stderr, "Division by zero! Exiting...");
	exit (EXIT_FAILURE);
    }

    if (c == NULL || AA_SIZE(c) == 0) {
	//c is zero
	*res_a = makePolynomial_AAZ (0, nvar);
	return;
    }
    
    // b is a monomial so we can do a simple divide
    if (AA_SIZE(b) == 1) {
	AltArrZ_t* res_r;
	divideBySingleTerm_AAZ (c, b, res_a, &res_r, nvar);
	//fprintf (stderr, "r := \n" );
	//printAAZ (res_r);
	freePolynomial_AAZ (res_r);
	return;
    }

    AAZElem_t* __restrict__ k = c->elems;
    AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
    AAZElem_t* __restrict__ b2Elem = b->elems + 1;
    register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
    // register int maxRSize = maxASize;
    register int i = 0;
    register int j = 0; 
    
    AltArrZ_t* a = makePolynomial_AAZ (maxASize, nvar);
    // AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
    AAZElem_t* __restrict__ curA = a->elems;
    // AAZElem_t* __restrict__ curR = r->elems;
    mpz_init(curA->coef);
    // mpz_init(curR->coef);
    
    register degrees_t beta = b->elems->degs;
    
    //init a with lt(c)/lt(b);
    DivTest_ptr divTest = getMonomialDivideTestFuncPtr(nvar);
    // while (k != lenK) { // && !( (*divTest)(k->degs, beta) && mpz_divisible_p(k->coef, b->elems->coef)) 
    // mpz_set(curR->coef, k->coef);
    // curR->degs = k->degs;
    // ++j;
    // if (j >= maxRSize) {
    // maxRSize <<= 1;
    //  r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
    //   curR = r->elems + j - 1; 	
    // }
    // ++(curR);
    // mpz_init(curR->coef);
    // ++k;
    // }
    
    if (k == lenK) {
	//no division to do at all!
	mpz_clear (curA->coef);
	// mpz_clear(curR->coef);
	
	AA_SIZE(a) = i;
	// AA_SIZE(r) = j;
	a->alloc = maxASize;
	// r->alloc = maxRSize;
	
	*res_a = a;
	// *res_r = r;
	return; 
    }

    subtractExponentVectors (k->degs, beta, curA->degs);
    //assuming b is monic!
    // mpz_set(curA->coef, k->coef);
    mpz_divexact (curA->coef, k->coef, b->elems->coef);
    ++k;
    
    //init multiplication between a (quotient) and b (divisor)
    ProductHeap_AAZ* h = prodheapCreate_AAZ (nvar);
    prodheapResize_AAZ (h, maxASize);
    h->lastB = AA_SIZE(b) - 1;
    prodheapInsert_AAZ (h, prodheapMakeChain_AAZ (0, 1, NULL), curA->degs + b2Elem->degs);
    ++i;
    ++curA;
    mpz_init (curA->coef);
    
    degrees_t* delta = prodheapPeek_AAZ (h); 
    register degrees_t eps;
    register cmpExp_t cmp;
    while (delta != NULL || k != lenK) {
	
	if (k == lenK) {
	    if (delta == NULL) {
		break;
	    }
	    cmp = 1;
	} else if (delta == NULL) {
	    cmp = -1;
	} else {
	    cmp = compareExponentVectors (*delta, k->degs);
	}
	
	if (cmp > 0) {
	    eps = *delta;
	    divisionGetNextTerm_AAZ (h, a->elems, b->elems, &(curA->coef));
	    if (mpz_sgn (curA->coef) == 0) {
		//in this case, the term with degree delta ended up 
		//having its coffeicient cancelled out (i.e. 0)
		delta = prodheapPeek_AAZ (h);
		continue;
	    } else {
		mpz_neg (curA->coef, curA->coef);
	    }
	} else if (cmp == 0) {
	    eps = *delta;
	    divisionGetNextTerm_AAZ (h, a->elems, b->elems, &(curA->coef));
	    if (mpz_sgn (curA->coef) == 0) {
		delta = prodheapPeek_AAZ (h);
		continue; //the chains cancelled themselves out since the peek
	    } else {
		mpz_sub (curA->coef, k->coef, curA->coef);
		++k;
		if (mpz_sgn (curA->coef) == 0) {
		    delta = prodheapPeek_AAZ (h);
		    continue;
		}
	    }
	} else {
	    eps = k->degs;
	    mpz_set (curA->coef,k->coef);
	    ++k;
	}
	
	if ((*divTest)(eps, beta) && mpz_divisible_p(curA->coef, b->elems->coef)) {
	    subtractExponentVectors(eps, beta, curA->degs);
	    mpz_divexact(curA->coef, curA->coef, b->elems->coef);
	    if (i+1 >= maxASize) {
		maxASize <<= 1;
		a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
		curA = a->elems + i;
		//prodheap maximum size should be equal to the size of a  	
		prodheapResize_AAZ(h, maxASize);
	    }
	    prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
	    ++i;
	    ++(curA);
	    mpz_init(curA->coef);
	}
	// else {
	    
	    //swap here so that curA becomes 0.
	    //mpz_swap(curR->coef, curA->coef);
	    //curR->degs = eps;
	    //++j;
	    /* if (j >= maxRSize) { */
	    /* 	maxRSize <<= 1; */
	    /* 	r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t)); */
	    /* 	curR = r->elems + j - 1; 	 */
	    /* } */
	    // ++(curR);
	    // mpz_init(curR->coef);
	// }
	
	delta = prodheapPeek_AAZ(h);
    }
    
    prodheapFree_AAZ(h);
    
    //clear since we always setup one past where we actually are.
    mpz_clear(curA->coef);
    //mpz_clear(curR->coef);
    
    AA_SIZE(a) = i;
    //AA_SIZE(r) = j;
    a->alloc = maxASize;
    //r->alloc = maxRSize;
    
    *res_a = a;
    //*res_r = r;
    return;
}

static void multiplyByInteger_AAZ_inp(AltArrZ_t* aa, const mpz_t z) {
	if (aa == NULL || aa->size == 0) {
		return;
	}

	for (int i = 0; i < aa->size; ++i) {
		mpz_mul(aa->elems[i].coef, aa->elems[i].coef, z);
	}
}

void univariatePseudoDivideBySingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
	if (b == NULL || b->size == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	if (c == NULL) {
		//c is zero
		*res_a = NULL;
		*res_r = NULL;
		return;
	}

	int nvar = 1;
	c = deepCopyPolynomial_AAZ(c);

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;
	register int j = 0; 

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	AAZElem_t* curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	int multSteps = 0;

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar)) {
			if (mpz_divisible_p(k->coef, bElems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, bElems->coef);
				multiplyByInteger_AAZ_inp(a, bElems->coef);
				++multSteps;
			}

			mpz_divexact(curA->coef, k->coef, bElems->coef);
			// mpz_set(curA->coef, k->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			mpz_set(curR->coef, k->coef);
			curR->degs = k->degs;
			++j;
			++(curR);
			mpz_init(curR->coef);
		}
		++k;
	}

	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i; 
	AA_SIZE(r) = j; 
	a->alloc = maxSize;
	r->alloc = maxSize;

	if (!lazy) {
		degrees_t d = c->elems->degs - b->elems->degs + 1 - multSteps;
		multSteps += d;

		mpz_t bPow;
		mpz_init(bPow);
		mpz_set(bPow, b->elems->coef);
		mpz_pow_ui(bPow, bPow, d);
		multiplyByInteger_AAZ_inp(a, bPow);
		multiplyByInteger_AAZ_inp(r, bPow);
		mpz_clear(bPow);

		// for (int j = 0; j < d; ++j) {
		// 	multiplyByInteger_AAZ_inp(a, b->elems->coef);
		// 	multiplyByInteger_AAZ_inp(r, b->elems->coef);
		// }
	}

	*res_a = a;
	*res_r = r;

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);
}

void univariatePseudoDividePolynomials_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, AltArrZ_t** res_r, int* e, int lazy) {
	if (b == NULL || AA_SIZE(b) == 0) {
		//division by zero
		fprintf(stderr, "Division by zero! Exiting...");
		exit(EXIT_FAILURE);
	}

	int nvar = 1;
	if (c == NULL || AA_SIZE(c) == 0) {
		//c is zero
		*res_a = makePolynomial_AAZ(0, nvar);
		*res_r = NULL;
		return;
	}

	if (isLessExponentVectors(c->elems->degs, b->elems->degs)) {
		*res_r = deepCopyPolynomial_AAZ(c);
		*res_a = NULL;
		if (e != NULL) {
			*e = 0;
		}
		return;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		univariatePseudoDivideBySingleTerm_AAZ(c, b, res_a, res_r, e, lazy);
		return;
	}

	c = deepCopyPolynomial_AAZ(c);

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register int maxASize = AA_SIZE(c) < 5 ? 5 : AA_SIZE(c);
	register int maxRSize = maxASize;
	register int i = 0;
	register int j = 0; 

	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AltArrZ_t* r = makePolynomial_AAZ(maxRSize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	AAZElem_t* __restrict__ curR = r->elems;
	mpz_init(curA->coef);
	mpz_init(curR->coef);

	register degrees_t beta = b->elems->degs;

	//init a with lt(c)/lt(b);
	int multSteps = 0;
	subtractExponentVectors(k->degs, beta, curA->degs);
	if (mpz_divisible_p(k->coef, b->elems->coef) == 0) {
		//set the quotient coefficient before updating dividend!
		mpz_set(curA->coef, k->coef);
		multiplyByInteger_AAZ_inp(c, b->elems->coef);
		++multSteps;
	} else {
		mpz_divexact(curA->coef, k->coef, b->elems->coef);
	}
	
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	a->size = 1;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h); 
	register degrees_t eps;
	register cmpExp_t cmp;
	while(delta != NULL || k != lenK) {
		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (compareExponentVectors(eps, beta) >= 0) {
			if (mpz_divisible_p(curA->coef, b->elems->coef) == 0) {
				multiplyByInteger_AAZ_inp(c, b->elems->coef);
				multiplyByInteger_AAZ_inp(a, b->elems->coef);
				++multSteps;
			
			} else {
				mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			}

			subtractExponentVectors(eps, beta, curA->degs);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i;
				//prodheap maximum size should be equal to the size of a  	
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			++i;
			++(curA);
			++(a->size);
			mpz_init(curA->coef);
		} else {
			//swap here so that curA becomes 0.
			mpz_swap(curR->coef, curA->coef);
			curR->degs = eps;
			++j;
			if (j >= maxRSize) {
				maxRSize <<= 1;
				r->elems = (AAZElem_t*) realloc(r->elems, maxRSize*sizeof(AAZElem_t));
				curR = r->elems + j - 1; 	
			}
			++(curR);
			mpz_init(curR->coef);
		}

		delta = prodheapPeek_AAZ(h);
	}

	prodheapFree_AAZ(h);

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);
	mpz_clear(curR->coef);

	AA_SIZE(a) = i;
	AA_SIZE(r) = j;
	a->alloc = maxASize;
	r->alloc = maxRSize;

	if (!lazy) {
		degrees_t d = c->elems->degs - b->elems->degs + 1 - multSteps;
		multSteps += d;

		mpz_t bPow;
		mpz_init(bPow);
		mpz_set(bPow, b->elems->coef);
		mpz_pow_ui(bPow, bPow, d);
		multiplyByInteger_AAZ_inp(a, bPow);
		multiplyByInteger_AAZ_inp(r, bPow);
		mpz_clear(bPow);

		// for (int j = 0; j < d; ++j) {
			// multiplyByInteger_AAZ_inp(a, b->elems->coef);
			// multiplyByInteger_AAZ_inp(r, b->elems->coef);
		// }
	}

	*res_a = a;
	*res_r = r;

	if (e != NULL) {
		*e = multSteps;
	}

	freePolynomial_AAZ(c);
}


int divideTestSingleTerm_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {
	if (b == NULL || b->size == 0) {
		return 0;
	}

	if (c == NULL || c->size == 0) {
		//c is zero
		*res_a = NULL;
		return 1;
	}

	AAZElem_t* k = c->elems;
	AAZElem_t* lenK = k + AA_SIZE(c);

	int maxSize = AA_SIZE(c) + 1;
	register int i = 0;

	AltArrZ_t* a = makePolynomial_AAZ(maxSize, nvar);
	AAZElem_t* bElems = b->elems;
	AAZElem_t* curA = a->elems;
	mpz_init(curA->coef);

	while (k != lenK) {
		if (monomialDivideTest(k->degs,bElems->degs,nvar) && mpz_divisible_p(k->coef, bElems->coef)) {
			mpz_divexact(curA->coef, k->coef, bElems->coef);
			subtractExponentVectors(k->degs, bElems->degs, curA->degs);
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {
			a->size = i;
			freePolynomial_AAZ(a);
			return 0;
		}
		++k;
	}

	mpz_clear(curA->coef);

	AA_SIZE(a) = i; 
	a->alloc = maxSize;

	*res_a = a;

	return 1;
}

int divideTest_AAZ(AltArrZ_t* c, AltArrZ_t* b, AltArrZ_t** res_a, int nvar) {
	
	if (b == NULL || AA_SIZE(b) == 0) {
		return 0;
	}

	if (c == NULL || AA_SIZE(c) == 0) {
		*res_a = makePolynomial_AAZ(0, nvar);
		return 1;
	}

	// b is a monomial so we can do a simple divide
	if (AA_SIZE(b) == 1) {
		return	divideTestSingleTerm_AAZ(c, b, res_a, nvar);
	}

	AAZElem_t* __restrict__ k = c->elems;
	AAZElem_t* __restrict__ lenK = k + AA_SIZE(c);
	AAZElem_t* __restrict__ b2Elem = b->elems + 1;
	register degrees_t beta = b->elems->degs;

	if (!monomialDivideTest(k->degs, beta, nvar) || !mpz_divisible_p(k->coef, b->elems->coef)) {
		return 0;
	}

	register int maxASize = AA_SIZE(c) + 1;
	register int i = 0;
	AltArrZ_t* a = makePolynomial_AAZ(maxASize, nvar);
	AAZElem_t* __restrict__ curA = a->elems;
	mpz_init(curA->coef);
	
	//init a with lt(c)/lt(b);
	subtractExponentVectors(k->degs, beta, curA->degs);
	mpz_divexact(curA->coef, k->coef, b->elems->coef);
	++k;

	//init multiplication between a (quotient) and b (divisor)
	ProductHeap_AAZ* h = prodheapCreate_AAZ(nvar);
	prodheapResize_AAZ(h, maxASize);
	h->lastB = AA_SIZE(b) - 1;
	prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(0, 1, NULL), curA->degs + b2Elem->degs);
	++i;
	++curA;
	mpz_init(curA->coef);

	degrees_t* delta = prodheapPeek_AAZ(h);
	register degrees_t eps;
	register cmpExp_t cmp;
	while(k != lenK || delta != NULL) {
		
		if (k == lenK) {
			if (delta == NULL) {
				break;
			}
			cmp = 1;
		} else if (delta == NULL) {
			cmp = -1;
		} else {
			cmp = compareExponentVectors(*delta, k->degs);
		}

		if (cmp > 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				//in this case, the term with degree delta ended up 
				//having its coffeicient cancelled out (i.e. 0)
				delta = prodheapPeek_AAZ(h);
				continue;
			} else {
				mpz_neg(curA->coef, curA->coef);
			}
		} else if (cmp == 0) {
			eps = *delta;
			divisionGetNextTerm_AAZ(h, a->elems, b->elems, &(curA->coef));
			if (mpz_sgn(curA->coef) == 0) {
				delta = prodheapPeek_AAZ(h);
				continue; //the chains cancelled themselves out since the peek
			} else {
				mpz_sub(curA->coef, k->coef, curA->coef);
				++k;
				if (mpz_sgn(curA->coef) == 0) {
					delta = prodheapPeek_AAZ(h);
					continue;
				}
			}
		} else {
			eps = k->degs;
			mpz_set(curA->coef,k->coef);
			++k;
		}

		if (monomialDivideTest(eps, beta, nvar) && mpz_divisible_p(curA->coef, b->elems->coef)) {
			subtractExponentVectors(eps, beta, curA->degs);
			mpz_divexact(curA->coef, curA->coef, b->elems->coef);
			if (i+1 >= maxASize) {
				maxASize <<= 1;
				a->elems = (AAZElem_t*) realloc(a->elems, maxASize*sizeof(AAZElem_t));
				curA = a->elems + i - 1;
				//prodheap maximum size should be equal to the size of a  	
				prodheapResize_AAZ(h, maxASize);
			}
			prodheapInsert_AAZ(h, prodheapMakeChain_AAZ(i, 1, NULL), curA->degs + b2Elem->degs);
			
			++i;
			++(curA);
			mpz_init(curA->coef);
		} else {

		    //divide test fails
		    
		    prodheapFree_AAZ(h);
		    a->size = i;
		    freePolynomial_AAZ(a);
		    return 0;
		}

		delta = prodheapPeek_AAZ(h);
	}

	//clear since we always setup one past where we actually are.
	mpz_clear(curA->coef);

	AA_SIZE(a) = i;
	a->alloc = maxASize;

	*res_a = a;
	return 1;
}




/*****************
 * Derivative / Integral
 *****************/

AltArrZ_t* derivative_AAZ(AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx < 0 || idx >= aa->nvar) {
		return NULL;
	}

	AltArrZ_t* ret = makePolynomial_AAZ(aa->size, aa->nvar);

	int* sizes = getExpOffsetArray(aa->nvar);
	unsigned long long int* masks = getExpMaskArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* __restrict__ elems = aa->elems;
	AAZElem_t* __restrict__ retElems = ret->elems;
	degrees_t deg;
	mpz_t mpzDeg;
	mpz_init(mpzDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		if (!(elems[i].degs & masks[idx])) {
			continue;
		}
		deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		if (deg < k) {
			continue;
		}

		retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpz_init(retElems[insertIdx].coef);
		mpz_set(retElems[insertIdx].coef, elems[i].coef);

		mpz_set_ui(mpzDeg, deg);
		for (int j = 0; j < k; ++j) {
			mpz_mul(retElems[insertIdx].coef, retElems[insertIdx].coef, mpzDeg);			
			--deg;	
			mpz_sub(mpzDeg, mpzDeg, mpzOne);
		}
		retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpz_clear(mpzDeg);
	mpz_clear(mpzOne);

	free(sizes);
	free(masks);

	ret->size = insertIdx;
	return ret;
}


/**
 * Integrate with respect to a variable that does not currently exist in aa.
 * It becomes the main variable (that is, to the left).
 *
 */
AltArr_t* integrateExpand_AAZ(AltArrZ_t* aa, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	unsigned long long int* maxs = getMaxExpArray((aa->nvar)+1);
	degrees_t degsK = k;
	if (degsK > maxs[0]) {
		fprintf(stderr, "SMQP exponent overflow in integration! At index 0, value: %d\n", k);
		exit(1);
	}
	free(maxs);


	AltArr_t* ret = deepCopyPolynomial_AAFromAAZ(aa);
	expandNumVarsLeft_AA(ret, (aa->nvar)+1);

	int size = aa->size;
	AAElem_t* elems = ret->elems;

	mpq_t kFact;
	mpq_init(kFact);
	mpz_fac_ui(mpq_numref(kFact), (unsigned long)k);

	int expOffset = getMVarExpOffset(ret->nvar);

	for (int i = 0; i < size; ++i) {
		mpq_div(elems[i].coef, elems[i].coef, kFact);
		elems[i].degs |= (degsK << expOffset);
	}

	mpq_clear(kFact);

	return ret;
}


AltArr_t* integral_AAZ(AltArrZ_t* aa, int idx, int k) {
	if (aa == NULL || aa->size == 0) {
		return NULL;
	}

	if (idx >= aa->nvar || idx < 0) {
		return integrateExpand_AAZ(aa, k);
	}

	AltArr_t* ret = makePolynomial_AA(aa->size, aa->nvar);

	int* sizes = getExpOffsetArray(aa->nvar);
	unsigned long long int* masks = getExpMaskArray(aa->nvar);
	unsigned long long int* maxs = getMaxExpArray(aa->nvar);

	int insertIdx = 0;
	int size = aa->size;
	AAZElem_t* __restrict__ elems = aa->elems;
	AAElem_t* __restrict__ retElems = ret->elems;
	degrees_t deg;
	mpq_t mpqDeg;
	mpq_init(mpqDeg);
	mpz_t mpzOne;
	mpz_init(mpzOne);
	mpz_set_si(mpzOne, 1l);
	for (int i = 0; i < size; ++i) {
		deg = GET_NTH_EXP(elems[i].degs, masks[idx], sizes[idx]);
		
		retElems[insertIdx].degs = (elems[i].degs & ~(masks[idx]));
		mpq_init(retElems[insertIdx].coef);
		mpz_set(mpq_numref(retElems[insertIdx].coef), elems[i].coef);

		mpq_set_ui(mpqDeg, deg+1, 1ul);
		for (int j = 0; j < k; ++j) {
			mpq_div(retElems[insertIdx].coef, retElems[insertIdx].coef, mpqDeg);			
			++deg;
			mpz_add(mpq_numref(mpqDeg), mpq_numref(mpqDeg), mpzOne);
		}

		if (deg > maxs[idx]) {
			fprintf(stderr, "SMZP exponent overflow in integration! Index: %d, value %llu\n", idx, deg);
			exit(1);
		}

		retElems[insertIdx].degs |= (deg << sizes[idx]);

		++insertIdx;
	}

	mpq_clear(mpqDeg);
	mpz_clear(mpzOne);

	free(sizes);
	free(masks);
	free(maxs);

	ret->size = insertIdx;
	return ret;
}



/*****************
 * Content, PrimitivePart, etc.
 *****************/

void integralContent_AAZ(AltArrZ_t* aa, mpz_t ret) {
	if (aa == NULL || aa->size <= 0) {
	    mpz_set_ui(ret, 1ul);
	    return;
	}

	AAZElem_t* elems = aa->elems;
	int size = aa->size;
	mpz_set_ui(ret, 1ul);
	
	mpz_t one;
	mpz_init(one);
	mpz_set_si(one, 1l);
	
	mpz_t cont; 
	mpz_init(cont);
	mpz_abs(cont, elems->coef);
	
	for (int i = 1; i < size; ++i) {
	    if (mpz_cmp(one, cont) == 0) {
		break;
	    }
	    mpz_gcd(cont, cont, elems[i].coef);
	}
	
	if (mpz_sgn(elems->coef) < 0 && mpz_sgn(cont) > 0) {
	    mpz_neg(cont, cont);
	}
	
	mpz_set(ret, cont);
	mpz_clear(cont);
	mpz_clear(one);
}

AltArrZ_t* integralContent_AAZ_polyOut (AltArrZ_t* aa){
    

}

AltArrZ_t* primitivePart_AAZ(AltArrZ_t* aa) {
    mpz_t content;
    mpz_init(content);
    
    integralContent_AAZ(aa, content);
    
    AltArrZ_t* res = deepCopyPolynomial_AAZ(aa);
    if (mpz_cmp_si(content, 1l) == 0) {
	return res;
    }
    
    AAZElem_t* elems = res->elems;
    int size = res->size;
    for (int i = 0; i < size; ++i) {
	mpz_divexact(elems[i].coef, elems[i].coef, content);
    }
    
    mpz_clear(content);
    return res;
}

AltArrZ_t* primitivePartAndContent_AAZ(AltArrZ_t* aa, mpz_t cont) {
	
	integralContent_AAZ(aa, cont);

	AltArrZ_t* res = deepCopyPolynomial_AAZ(aa);
	if (mpz_cmp_si(cont, 1l) == 0) {
		return res;
	}

	AAZElem_t* elems = res->elems;
	int size = res->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, cont);
	}

	return res;
}

void primitivePart_AAZ_inp(AltArrZ_t* aa) {
	if (aa == NULL || aa->size <= 0) {
		return;
	}

	mpz_t content;
	mpz_init(content);

	integralContent_AAZ(aa, content);
	if (mpz_cmp_si(content, 1l) == 0) {
		mpz_clear(content);
		return;
	}

	AAZElem_t* elems = aa->elems;
	int size = aa->size;
	for (int i = 0; i < size; ++i) {
		mpz_divexact(elems[i].coef, elems[i].coef, content);
	}

	mpz_clear(content);
}

AltArrZ_t* univariateGCD_AAZ(AltArrZ_t* a, AltArrZ_t* b) {
	if (a->nvar != 1 || b->nvar != 1) {
		fprintf(stderr, "SMQP ERROR: Calling univariate GCD on multivariate arrays\n");
		exit(1);
	}


	if (a == NULL || a->size == 0) {
		return deepCopyPolynomial_AAZ(b);
	}
	if (b == NULL || b->size == 0){
		return deepCopyPolynomial_AAZ(a);
	}

	AltArrZ_t* r0 = NULL;
	AltArrZ_t* r1 = NULL;
	AltArrZ_t* r2 = NULL;

	mpz_t c0;
	mpz_t c1;
	mpz_init(c0);
	mpz_init(c1);

	if (isGreaterExponentVectors(a->elems->degs, b->elems->degs)) {
		r0 = primitivePartAndContent_AAZ(a, c0);
		r1 = primitivePartAndContent_AAZ(b, c1);
		// r0 = deepCopyPolynomial_AA(a);
		// r1 = deepCopyPolynomial_AA(b);
	} else {
		r0 = primitivePartAndContent_AAZ(b, c0);
		r1 = primitivePartAndContent_AAZ(a, c1);
		// r0 = deepCopyPolynomial_AA(b);
		// r1 = deepCopyPolynomial_AA(a);
	}

	AltArrZ_t* quo = NULL;
	while (r1 != NULL && r1->size > 0) {
		univariatePseudoDividePolynomials_AAZ(r0, r1, &quo, &r2, NULL, 1);
		freePolynomial_AAZ(quo);
		quo = NULL;

		freePolynomial_AAZ(r0);
		r0 = r1;
		r1 = r2;
		primitivePart_AAZ_inp(r1);
		r2 = NULL;
 	}

	freePolynomial_AAZ(r1);
	freePolynomial_AAZ(r2);
	if (r0 != NULL && r0->size > 0 && mpz_sgn(r0->elems->coef) < 0) {
		negatePolynomial_AAZ(r0);
	}

	mpz_clear(c0);
	mpz_clear(c1);

	return r0;

}

AltArrZ_t* commonFactor_AAZ(AltArrZ_t* a, AltArrZ_t** factored) {
	if (a == NULL || a->size == 0) {
		return NULL;
	}

	AltArrZ_t* ret = makePolynomial_AAZ(1, a->nvar);
	mpz_init(ret->elems->coef);
	mpz_set_ui(ret->elems->coef, 1ul);

	if (a->nvar == 1) {
		ret->elems->degs = a->elems[a->size-1].degs;
	} else {
		degrees_t min = a->elems[a->size-1].degs;
		degrees_t* masks = getExpMaskArray(a->nvar);
		for (int i = a->size-2; i >= 0; --i) {
			for (int j = 0; j < a->nvar; ++j) {
				if ((a->elems[i].degs & masks[j]) < (min & masks[j]) ) {
					min = min & (~masks[j]); //zero out j;
					min |= (a->elems[i].degs & masks[j]);					
				}
			}
			if (isZeroExponentVector(min)) {
				break;
			}
		} 
		free(masks);
		ret->elems->degs = min; 
	}
	ret->size = 1;
	
	if (factored != NULL) {
		AltArrZ_t* factRet = deepCopyPolynomial_AAZ(a);
		if (! isZeroExponentVector(ret->elems->degs)) {
			for (int i = 0; i < a->size; ++i) {
				factRet->elems[i].degs -= ret->elems->degs;
			}
		}
		*factored = factRet;
	}

	return ret;
}



/*****************
 * Evaluation
 *****************/

void univarEvaluate_AAZ(AltArrZ_t* aa, const mpz_t point, mpz_t res) {
	if (aa->size == 0) {
		mpz_set_si(res, 0l);
		return;
	}

	AAZElem_t* elems = aa->elems;
	register int size = aa->size;
	mpz_set(res, elems->coef);
	degrees_t prevDeg = elems->degs;
	degrees_t nextDeg;
	for (int i = 1; i < size; ++i) {
		nextDeg = elems[i].degs;
		for (degrees_t j = prevDeg; j > nextDeg; --j) {
			mpz_mul(res, res, point);
		}
		mpz_add(res, res, elems[i].coef);
		prevDeg = nextDeg;
	}
	for (degrees_t j = prevDeg; j > 0; --j) {
		mpz_mul(res, res, point);
	}
}



/****************
* Multi-Divisor Division
*****************/

void heapMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
  if (s < 0)
  {
    fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
    return;
  }

  if (s == 0)
  {
    *r = deepCopyPolynomial_AAZ (f);
    return;
  }

  if (s == 1)  // simple multivariate polynomial division
  {
	  if (G[0] != NULL && G[0]->size != 0){
		  dividePolynomials_AAZ (f, *G, Q, r, nvar);
	  } else {
		  *r = deepCopyPolynomial_AAZ (f);
		  Q[0] = NULL;
	  }
	  return;
  }

  int i; // i= 0 ... (s-1)
  AltArrZ_t* h = deepCopyPolynomial_AAZ (f);  // this algorithm runs until h != NULL
  AltArrZ_t* rem = NULL;
  AltArrZ_t* lt = NULL;
  AltArrZ_t* tmp_q = NULL;
  AltArrZ_t* tmp_r = NULL;

  while (h->size != 0 && h != NULL)
  {
    i = 0; // in each step i should start from the first divisor
    while (i < s)
    {
      if ( h->size != 0 && h!= NULL && G[i] != NULL && G[i]->size != 0  && monomialDivideTest (h->elems[0].degs, G[i]->elems[0].degs, nvar))
      {
      	tmp_q = NULL;
      	tmp_r = NULL;
        dividePolynomials_AAZ (h, G[i], &tmp_q, &tmp_r, nvar);
        if (Q[i] != NULL)
	        Q[i] = addPolynomials_AAZ (Q[i], tmp_q, nvar);
	    else
	    	Q[i] = tmp_q;
        h = tmp_r;
        i = 0;
      } else
        i = i + 1;
    }
    if (h->size != 0 && h != NULL){
      lt = leadingTerm_AAZ (h, nvar);
    }
    else
      lt = NULL;

    if (rem == NULL)
      rem = lt;
    else
    {
      if (lt != NULL) {
        rem = addPolynomials_AAZ_inp (rem, lt, nvar);
      }
    }
    if (h->size != 0 && h!= NULL)
        h = subPolynomials_AAZ_inp (h, lt, nvar);
  }

  *r = rem;
  return;
}


void triangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
	if (s < 0)
	{
		fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
		return;
	}

	if (s == 0)
	{
		*r = deepCopyPolynomial_AAZ (f);
		return;
	}

	if (s == 1)  // simple multivariate polynomial division
	{
		  if (G[0] != NULL && G[0]->size != 0){
			  dividePolynomials_AAZ (f, *G, Q, r, nvar);
		  } else {
			  *r = deepCopyPolynomial_AAZ (f);
			  Q[0] = NULL;
		  }
		  return;
	}

	int index = 0;
	int i = 0;
	int* orderList = recursiveLoop(s);
	int orderSize = (1 << s) - 1;
	AltArrZ_t* h =deepCopyPolynomial_AAZ (f);
	AltArrZ_t* rem = NULL;
	AltArrZ_t* lt = NULL;
	AltArrZ_t* tmp_q = NULL;
	AltArrZ_t* tmp_r = NULL;

	while (h->size != 0 && h != NULL)
	{
		for (i = 0; i < orderSize; i++)
		{
			index = orderList[i];
			if (h->size != 0 && h != NULL && G[index] != NULL  && G[index]->size != 0 && monomialDivideTest(h->elems[0].degs , G[index]->elems[0].degs , nvar))
			{
				tmp_q = NULL;
				tmp_r = NULL;

				dividePolynomials_AAZ (h, G[index], &tmp_q, &tmp_r, nvar);
				if (Q[index] != NULL) {
					Q[index] = addPolynomials_AAZ_inp(Q[index], tmp_q, nvar);
				}
				else {
					Q[index] = tmp_q;
				}
				freePolynomial_AAZ (h);
				h = tmp_r;
			}
		}
		if (h->size != 0 && h != NULL)
			lt = leadingTerm_AAZ (h, nvar);
		else
			lt = NULL;

		if (rem == NULL)
			rem = lt;
		else {
			if (lt != NULL){
				rem = addPolynomials_AAZ_inp (rem, lt, nvar);
			}
		}
		if (h != NULL && h->size != 0)
			h = subPolynomials_AAZ_inp (h, lt, nvar);
	}

	*r = rem;
	free(orderList);
	return;
}

// TODO:
void primitiveFactorTriangularSetMDD_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar)
{
  if (s < 0)
  {
    fprintf(stderr, "BPAS Error: the number of divisor set is out of range!");
    return;
  }

  if (s == 0)
  {
    *r = deepCopyPolynomial_AAZ (f);
    return;
  }

  if (s == 1)  // simple multivariate polynomial division
  {
	  if (G[0] != NULL && G[0]->size != 0){
		  dividePolynomials_AAZ (f, *G, Q, r, nvar);
	  } else {
		  *r = deepCopyPolynomial_AAZ (f);
		  Q[0] = NULL;
	  }
    return;
  }

  int index = 0;
  int i = 0;
  int* orderList = recursiveLoop(s);
  int orderSize = (1 << s) - 1;
  AltArrZ_t* h =deepCopyPolynomial_AAZ (f);
  AltArrZ_t* rem = NULL;
  AltArrZ_t* lt = NULL;
  AltArrZ_t* tmp_q = NULL;
  AltArrZ_t* tmp_r = NULL;

  while (h->size != 0 && h != NULL)
  {
    for (i = 0; i < orderSize; i++)
    {
      index = orderList[i];
      if (h->size != 0 && h != NULL && G[index] != NULL && G[index]->size != 0  && monomialDivideTest(h->elems[0].degs , G[index]->elems[0].degs , nvar))
      {
      	tmp_q = NULL;
      	tmp_r = NULL;
        dividePolynomials_AAZ (h, G[index], &tmp_q, &tmp_r, nvar);
        if (Q[index] != NULL)
	        Q[index] = addPolynomials_AAZ_inp(Q[index], tmp_q, nvar);
	    else
	    	Q[index] = tmp_q;
        h = tmp_r;
      }
    }
    if (h->size != 0 && h != NULL)
      lt = leadingTerm_AAZ (h, nvar);
    else
      lt = NULL;

    if (rem == NULL)
      rem = lt;
    else {
      if (lt != NULL){
          rem = addPolynomials_AAZ_inp (rem, lt, nvar);
        }
    }
    if (h->size != 0 && h != NULL)
        h = subPolynomials_AAZ_inp (h, lt, nvar);
  }

  *r = rem;
  free(orderList);
  return;
}


void multiDivisorDivision_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t** r, int s, int nvar, int type)
{
	if (type == 0){
//		printf("Running HeapMDD ...\n");
		heapMDD_AAZ (f, G, Q, r, s, nvar);
		return;
	} else if (type == 1){
//		printf("Running triangularSetMDD ...\n");
		triangularSetMDD_AAZ (f, G, Q, r, s, nvar);
		return;
	} else{
//		printf("Running primitiveFactorTriangularSetMDD ...\n");
		primitiveFactorTriangularSetMDD_AAZ (f, G, Q, r, s, nvar);
	}
	return;
}


int multiDivisorDivisionVerification_AAZ (AltArrZ_t* f, AltArrZ_t** G, AltArrZ_t** Q, AltArrZ_t* r, AltArrZ_t* hPow, int nSet, int tnvar)
{
    int nullity = 0;
    for (int i = 0; i < nSet; ++i){
        if (G[i] == NULL || G[i]->size == 0){
            fprintf(stderr, "divisor[%d] is NULL!\n", i);
            nullity++;
        }
        if (Q[i] == NULL || Q[i]->size == 0){
            fprintf(stderr, "quotient[%d] is NULL!\n", i);
            nullity++;
        }
    }
    if (r == NULL || r->size == 0){
        fprintf(stderr, "remainder is NULL!\n");
        nullity++;
    }
    if (hPow == NULL || hPow->size == 0){
        fprintf(stderr, "hPow is NULL!\n");
        nullity++;
    }
    if (f == NULL || f->size == 0){
        fprintf(stderr, "dividend is NULL!\n");
        nullity++;
        exit(1);
    }
    

    AltArrZ_t* sum = NULL;
    for (int i = 0; i < nSet; ++i){
        if (G[i] != NULL && G[i]->size != 0 && Q[i] != NULL && Q[i]->size != 0){
            if (sum == NULL || sum->size == 0){
                sum = multiplyPolynomials_AAZ (G[i], Q[i], tnvar);
            } else {
            	sum = addPolynomials_AAZ_inp (sum , multiplyPolynomials_AAZ (G[i], Q[i], tnvar), tnvar);
            }
        }
    }
    if (r != NULL && r->size != 0){
        sum = addPolynomials_AAZ_inp (sum, r, tnvar);
    }

    AltArrZ_t* hf = NULL;
    if (hPow != NULL && hPow->size != 0){
         hf = multiplyPolynomials_AAZ (f, hPow, tnvar); 
    } else {
        hf = deepCopyPolynomial_AAZ (f);
    }

    AltArrZ_t* res = subPolynomials_AAZ (sum, hf, tnvar);
    if (res == NULL || res->size == 0){
        return 1;
    } else {
        fprintf(stderr, "No. Zero-Poly(s): %d\n", nullity);
        printf("The difference: \n");
        printAAZ (res);
        printf("\n");
        return 0;
    }
}

