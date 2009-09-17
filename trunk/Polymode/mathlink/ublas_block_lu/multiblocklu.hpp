template<class T>
class TridiagonalBlockLU
{
private:
	bool overwrite;
	T shift;
	ublas::matrix<std::size_t> pivots_m;
	block_array_wrapped<T> LU;
	block_array_wrapped<T> LUupdate;

public:
	//Perform Block LU
	bool lu(block_array_wrapped<T>& A, T shift);
	//Link to other class's LU decompostion
	bool link_lu(tridiagonal_block_lu<T>& tblu)
	{
	};

	bool lu_update(block_array_wrapped<T>& A, int nrowsup);

	//Back subsistution
	template<class M, class S>
	void solve(ublas::matrix<T,M,S>& x);

	//Back subsistution of transpose problem
	template<class M, class S>
	void solve_transpose(ublas::matrix<T,M,S>& x);
};

//Example:
	lua = TridiagonalBlockLu();
	lua.lu(A, shift);

//Create other linked LU objects
	lub = TridiagonalBlockLu(lua);
	luc = TridiagonalBlockLu(lub);

//But if lua.lu is updated, this information is lost!

//Or.. new class that uses the linked object
	template<class T>
	class LinkedBlockLU
	{
	private:
		TridiagonalBlockLU* linklu;
	public:
		//Create link to LU class
		LinkedBlockLU(tridiagonal_block_lu<T>& lu_in)
		{
			linklu = &lu_in;
		};
	...
	}

	lua = TridiagonalBlockLu();
	lub = LinkedBlockLu(lua);
	...
	lua.lu(A, shift);

//But now is lua is deleted the linked classes don't know!



