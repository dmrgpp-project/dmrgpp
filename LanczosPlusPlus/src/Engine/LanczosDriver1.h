#ifndef LANCZOSDRIVER_1_H
#define LANCZOSDRIVER_1_H
#include "LabeledOperator.h"
#include "LanczosDriver.h"
#include "LanczosGlobals.h"
#include <PsimagLite/Io/IoSimple.h>

template <typename ModelType> SizeType maxOrbitals(const ModelType& model)
{
	SizeType res = 0;
	for (SizeType i = 0; i < model.geometry().numberOfSites(); i++) {
		if (res < model.orbitals(i))
			res = model.orbitals(i);
	}
	return res;
}

template <typename EngineType>
void extendedStatic(PsimagLite::String                   manypoint,
                    const EngineType&                    engine,
                    const typename EngineType::PairType& braAndKet)
{
	typedef typename EngineType::VectorSizeType  VectorSizeType;
	PsimagLite::Vector<PsimagLite::String>::Type str;
	PsimagLite::split(str, manypoint, ";");

	VectorSizeType                                             sites;
	VectorSizeType                                             spins;
	VectorSizeType                                             orbs;
	PsimagLite::Vector<LanczosPlusPlus::LabeledOperator>::Type whats;
	for (SizeType i = 0; i < str.size(); ++i) {
		PsimagLite::Vector<PsimagLite::String>::Type str2;
		PsimagLite::split(str2, str[i], "?");
		if (str2.size() < 3)
			throw PsimagLite::RuntimeError("-M option malformed\n");
		whats.push_back(LanczosPlusPlus::LabeledOperator(str2[0]));
		sites.push_back(atoi(str2[1].c_str()));
		spins.push_back(atoi(str2[2].c_str()));
		if (str2.size() == 4)
			orbs.push_back(atoi(str2[3].c_str()));
		else
			orbs.push_back(0);
	}

	std::cout << "<gs|" << manypoint << "|gs>=";
	std::cout << engine.manyPoint(sites, whats, spins, orbs, braAndKet);
	std::cout << "\n";
}

template <typename ModelType,
          typename SpecialSymmetryType,
          template <typename, typename> class InternalProductTemplate>
void mainLoop3(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType         GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef LanczosPlusPlus::Engine<ModelType, InternalProductTemplate, SpecialSymmetryType>
	    EngineType;
	using ContinuedFractionCollectionType = PsimagLite::ContinuedFractionCollection<RealType>;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	const GeometryType& geometry = model.geometry();
	EngineType          engine(model, io);

	//! get the g.s.:
	RealType eg = engine.energies(0);
	std::cout.precision(8);
	std::cout << "Energy=" << eg << "\n";
	PsimagLite::String filename = PsimagLite::basenameOf(io.filename());

	const SizeType nmeas = lanczosOptions.measure.size();
	for (SizeType i = 0; i < nmeas; ++i) {
		VectorStringType tokens;
		PsimagLite::split(tokens, lanczosOptions.measure[i], ",");
		const SizeType ntokens = tokens.size();
		for (SizeType j = 0; j < ntokens; ++j) {
			VectorStringType bra_op_ket;
			PsimagLite::split(bra_op_ket, tokens[j], "|");
			engine.measure(bra_op_ket);
		}
	}

	bool needs_dos = false;
	try {
		int tmp = 0;
		io.readline(tmp, "ComputeDensityOfStates=");
		needs_dos = (tmp > 0);
	} catch (std::exception&) { }

	typedef std::pair<SizeType, SizeType>  PairSizeType;
	PsimagLite::Vector<PairSizeType>::Type pair_of_sites;
	const SizeType                         n = geometry.numberOfSites();

	if (needs_dos) {
		lanczosOptions.gf.push_back(LanczosPlusPlus::LabeledOperator("c"));
		for (SizeType i = 0; i < n; ++i)
			pair_of_sites.push_back(PairSizeType(i, i));
	}

	try {
		io.read(lanczosOptions.sites, "TSPSites");

		if (lanczosOptions.sites.size() == 0)
			err("TSPSites must have at least one site\n");

		if (lanczosOptions.sites.size() == 1)
			lanczosOptions.sites.push_back(lanczosOptions.sites[0]);

		pair_of_sites.push_back(
		    PairSizeType(lanczosOptions.sites[0], lanczosOptions.sites[1]));
	} catch (std::exception&) { }

	bool     has_center  = false;
	SizeType center_site = 0;
	try {
		io.readline(center_site, "TSPCenter=");
		std::cout << "TSPCenter=" << center_site << "\n";

		for (SizeType i = 0; i < n; ++i)
			pair_of_sites.push_back(PairSizeType(center_site, i));
		has_center = true;
	} catch (std::exception&) { }

	bool do_all_pairs = false;
	try {
		int tmp = 0;
		io.readline(tmp, "DoAllPairs=");
		do_all_pairs = (tmp > 0);
	} catch (std::exception&) { }

	if (do_all_pairs && has_center)
		err("You cannot have both TSPCenter and DoAllPairs\n");

	if (do_all_pairs) {
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = 0; j < n; ++j)
				pair_of_sites.push_back(PairSizeType(i, j));
	}

	for (SizeType gfi = 0; gfi < lanczosOptions.gf.size(); ++gfi) {
		SizeType       counter   = 0;
		const SizeType n_indices = pairOfSites.size();
		for (SizeType s_index = 0; s_index < n_indices; ++s_index) {
			const SizeType site0 = pairOfSites[s_index].first;
			const SizeType site1 = pairOfSites[s_index].second;

			std::cout << "#gf(i=" << site0 << ", j=" << site1 << ")\n";

			typename EngineType::VectorStringType vstr;
			PsimagLite::IoSimple::Out io_out(filename + ttos(counter) + ".comb");

			io_out.write(site0, "Site0");
			io_out.write(site1, "Site1");

			if (has_center)
				io_out.write(center_site, "TSPCenter");

			ContinuedFractionCollectionType cf_collection(PsimagLite::FreqEnum::REAL);
			SizeType                        norbitals = maxOrbitals(model);
			for (SizeType orb1 = 0; orb1 < norbitals; orb1++) {
				for (SizeType orb2 = orb1; orb2 < norbitals; orb2++) {
					engine.spectralFunction(
					    cfCollection,
					    vstr,
					    lanczosOptions.gf[gfi],
					    site0,
					    site1,
					    lanczosOptions.spins,
					    std::pair<SizeType, SizeType>(orb1, orb2));
				}
			}

			io_out << "#INDEXTOCF ";
			for (SizeType i = 0; i < vstr.size(); ++i)
				io_out << vstr[i] << " ";
			io_out << "\n";
			cf_collection.write(io_out);
			std::cerr << "LanczosDriver1.h: Written to " << io_out.filename() << "\n";
			++counter;
		}
	}

	for (SizeType cicji = 0; cicji < lanczosOptions.cicj.size(); cicji++) {
		SizeType                              total = geometry.numberOfSites();
		PsimagLite::Matrix<ComplexOrRealType> cicj_matrix(total, total);
		SizeType                              norbitals = maxOrbitals(model);
		for (SizeType orb1 = 0; orb1 < norbitals; orb1++) {
			for (SizeType orb2 = 0; orb2 < norbitals; orb2++) {
				engine.twoPoint(cicj_matrix,
				                lanczosOptions.cicj[cicji],
				                lanczosOptions.spins,
				                std::pair<SizeType, SizeType>(orb1, orb2),
				                std::pair<SizeType, SizeType>(0, 0));
				std::cout << cicj_matrix;
			}
		}
	}

	if (lanczosOptions.split >= 0) {
		LanczosPlusPlus::ReducedDensityMatrix<ModelType> reduced_density_matrix(
		    model, engine.eigenvector(0), lanczosOptions.split);
		reduced_density_matrix.printAll(std::cout);
	}

	if (lanczosOptions.extendedStatic != "") {
		PsimagLite::Vector<PsimagLite::String>::Type str;
		PsimagLite::split(str, lanczosOptions.extendedStatic, ",");
		for (SizeType i = 0; i < str.size(); ++i)
			extendedStatic(str[i], engine, typename EngineType::PairType(0, 0));
	}
}

template <typename ModelType, typename SpecialSymmetryType>
void mainLoop2(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	PsimagLite::String tmp;
	io.readline(tmp, "SolverOptions=");
	bool onthefly = (tmp.find("InternalProductOnTheFly") != PsimagLite::String::npos);

	if (onthefly) {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductOnTheFly>(
		    model, io, lanczosOptions);
	} else {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductStored>(
		    model, io, lanczosOptions);
	}
}

#endif // LANCZOSDRIVER_1_H
