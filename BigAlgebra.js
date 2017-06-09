class BigAlgebra extends BigArith{
	constructor(){
		super();
		this.ObjName = "BigAlgebra";
	}
	
	/**	Name of object 
	*	@return {String} - Name of the object
	*/
	get name(){
		return this.ObjName;
	}
	
	/**	Solve a simultaneous NxN equation
	*	function simultaneousEq
	*	@param eq - {Array} - A Nx[N+1] array. The +1 element representing the value after "="
	*	@param toFixed - {Number} - Optional. Fixed the result to @param decimal places
	*	@returns {Array} - An array of N elements representing each unknown 
	*/
	simultaneousEq(eq, toFixed){
		let eqlen = eq.length;
		//Check if equation is NxN
		for(let i = 0; i < eqlen; i++){
			if(eq[i].length != eqlen+1) throw new Error("Equation should be NxN");
		}
		
		let mat = [];
		//Get determinant of main equation
		for(let i = 0; i < eqlen; i++){
			mat.push(eq[i].slice(0, eq[i].length-1));
		}
		let D_main = this.matrixDet(mat);
		
		//Get determinant of each unknown
		let D_unknown = [];
		for(let i = 0; i < eqlen; i++){
			mat = [];
			for(let j = 0; j < eqlen; j++){
				let x = [].concat(eq[j].slice(0, i), eq[j][eq.length], eq[j].slice(i+1, eq[i].length-1));
				mat.push(x);
			}
			D_unknown.push(this.matrixDet(mat));
		}
		
		//Get results of unknowns
		let result = [];
		for(let i = 0; i < eqlen; i++){
			if(typeof toFixed != "undefined")
				result.push(new BigArith(D_unknown[i]/D_main).toFixed(toFixed).valueOf());
			else	
				result.push(D_unknown[i]/D_main);
		}
		return result;
	}
	
	/**	Find the determinant of a NxN matric
	*	function matrixDet
	*	@param {Array} - A NxN array.
	*	@returns {Number} - The determinant of @param 
	*/
	matrixDet(m){
		let mlen = m.length;
		//Check if matrix is NxN
		for(let i = 0; i < mlen; i++){
			if(m[i].length != mlen) throw new Error("Matrix should be NxN");
		}
		
		//Find determinant
		let det = 0;
		if(mlen > 2){
			let flag = true;
			for(let i = 0; i < mlen; i++){
				let mi = [];
				for(var t = 0; t < mlen; t++){
					let sub = [];
					for(var tt = 0; tt < mlen; tt++){
						if(t == 0 || tt == i) 
							continue;
						sub.push(m[t][tt]);
					}
					if(sub.length > 0) mi.push(sub);
				}
				if(flag){
					det += (m[0][i]*this.matrixDet(mi));
					flag = false;
				}
				else{
					det -= (m[0][i]*this.matrixDet(mi))
					flag = true;
				}
			}
		}
		else if(mlen == 2) det = (m[0][0]*m[1][1])-(m[0][1]*m[1][0]);
		else det = m[0][0];
		return det;
	}
}