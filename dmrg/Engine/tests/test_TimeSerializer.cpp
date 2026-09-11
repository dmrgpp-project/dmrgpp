#include <PsimagLite/PsimagLite.h>

#include "TimeSerializer.h"
#include <PsimagLite/Io/IoSelector.h>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <filesystem>
#include <string>
#include <vector>

namespace {

class TemporaryFile {
public:

	TemporaryFile()
	{
		const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
		path_            = std::filesystem::temp_directory_path()
		    / ("dmrgpp-time-serializer-" + std::to_string(stamp) + ".hdf5");
	}

	~TemporaryFile() { std::filesystem::remove(path_); }

	const std::filesystem::path& path() const { return path_; }

private:

	std::filesystem::path path_;
};

class TestVector {
public:

	using value_type = double;

	TestVector()
	    : value_(0)
	{ }

	explicit TestVector(double value)
	    : value_(value)
	{ }

	void write(PsimagLite::IoSelector::Out& io, PsimagLite::String label) const
	{
		io.write(value_, label);
	}

	void read(PsimagLite::IoSelector::In& io, PsimagLite::String label)
	{
		io.read(value_, label);
	}

	double value() const { return value_; }

private:

	double value_;
};

class TestAoe {
public:

	using VectorStageEnumType = std::vector<Dmrg::StageEnum>;

	TestAoe()
	    : vectors_({ TestVector(1.5), TestVector(2.5), TestVector(3.5) })
	    , stages_({ Dmrg::StageEnum::OPERATOR,
	                Dmrg::StageEnum::WFT_NOADVANCE,
	                Dmrg::StageEnum::WFT_ADVANCE })
	{ }

	SizeType tvs() const { return vectors_.size(); }

	const TestVector& targetVectors(SizeType index) const { return vectors_[index]; }

	const VectorStageEnumType& stages() const { return stages_; }

private:

	std::vector<TestVector> vectors_;
	VectorStageEnumType     stages_;
};

using TimeSerializerType = Dmrg::TimeSerializer<TestVector>;

void writeLegacyCheckpoint(const std::filesystem::path& path, const TestAoe& aoe)
{
	PsimagLite::IoSelector::Out io(path.string(), PsimagLite::IoSelector::ACC_TRUNC);
	io.createGroup("FinalPsi");
	TimeSerializerType serializer(4, 1.25, 2, aoe, "Expression");
	serializer.write(io, "FinalPsi");
	io.close();
}

void addPvectorTimes(const std::filesystem::path& path,
                     int                          version,
                     const std::vector<double>*   values)
{
	PsimagLite::IoSelector::Out io(path.string(), PsimagLite::IoSelector::ACC_RDW);
	const PsimagLite::String    prefix("FinalPsi/TimeSerializer/PvectorTimes");
	io.createGroup(prefix);
	io.write(version, prefix + "/Version");
	if (values)
		io.write(*values, prefix + "/Values");
	io.close();
}

} // namespace

TEST_CASE("TimeSerializer round-trips P-vector times", "[TimeSerializer]")
{
	TemporaryFile             temporaryFile;
	TestAoe                   aoe;
	const std::vector<double> pvectorTimes({ 0.0, 0.625, 0.7 });

	{
		PsimagLite::IoSelector::Out io(temporaryFile.path().string(),
		                               PsimagLite::IoSelector::ACC_TRUNC);
		io.createGroup("FinalPsi");
		TimeSerializerType serializer(4, 1.25, 2, aoe, "Expression", pvectorTimes);
		serializer.write(io, "FinalPsi");
		io.close();
	}

	PsimagLite::IoSelector::In io(temporaryFile.path().string());
	TimeSerializerType         serializer(io, "FinalPsi");

	REQUIRE(serializer.hasPvectorTimes());
	REQUIRE(serializer.pvectorTimes().size() == pvectorTimes.size());
	for (SizeType i = 0; i < pvectorTimes.size(); ++i) {
		CHECK(serializer.pvectorTime(i) == pvectorTimes[i]);
		CHECK(serializer.vector(i).value() == aoe.targetVectors(i).value());
	}
	CHECK(serializer.currentTimeStep() == 4);
	CHECK(serializer.time() == 1.25);
	CHECK(serializer.site() == 2);
	CHECK(serializer.name() == "Expression");
	REQUIRE(serializer.stages().size() == aoe.stages().size());
	for (SizeType i = 0; i < aoe.stages().size(); ++i)
		CHECK(serializer.stages()[i] == aoe.stages()[i]);
}

TEST_CASE("TimeSerializer reads legacy checkpoints without P-vector times", "[TimeSerializer]")
{
	TemporaryFile temporaryFile;
	TestAoe       aoe;
	writeLegacyCheckpoint(temporaryFile.path(), aoe);

	PsimagLite::IoSelector::In io(temporaryFile.path().string());
	TimeSerializerType         serializer(io, "FinalPsi");

	CHECK_FALSE(serializer.hasPvectorTimes());
	CHECK(serializer.pvectorTimes().empty());
	CHECK_THROWS(serializer.pvectorTime(0));
	CHECK(serializer.time() == 1.25);
	CHECK(serializer.numberOfVectors() == aoe.tvs());
}

TEST_CASE("TimeSerializer rejects unsupported P-vector-time versions", "[TimeSerializer]")
{
	TemporaryFile             temporaryFile;
	TestAoe                   aoe;
	const std::vector<double> values({ 0.0, 0.625, 0.7 });
	writeLegacyCheckpoint(temporaryFile.path(), aoe);
	addPvectorTimes(temporaryFile.path(), 2, &values);

	PsimagLite::IoSelector::In io(temporaryFile.path().string());
	CHECK_THROWS(TimeSerializerType(io, "FinalPsi"));
}

TEST_CASE("TimeSerializer rejects missing P-vector-time values", "[TimeSerializer]")
{
	TemporaryFile temporaryFile;
	TestAoe       aoe;
	writeLegacyCheckpoint(temporaryFile.path(), aoe);
	addPvectorTimes(temporaryFile.path(), 1, nullptr);

	PsimagLite::IoSelector::In io(temporaryFile.path().string());
	CHECK_THROWS(TimeSerializerType(io, "FinalPsi"));
}

TEST_CASE("TimeSerializer rejects P-vector-time count mismatch", "[TimeSerializer]")
{
	TemporaryFile             temporaryFile;
	TestAoe                   aoe;
	const std::vector<double> values({ 0.0, 0.625 });
	writeLegacyCheckpoint(temporaryFile.path(), aoe);
	addPvectorTimes(temporaryFile.path(), 1, &values);

	PsimagLite::IoSelector::In io(temporaryFile.path().string());
	CHECK_THROWS(TimeSerializerType(io, "FinalPsi"));
}

TEST_CASE("TimeSerializer rejects P-vector-time count mismatch before writing", "[TimeSerializer]")
{
	TestAoe                   aoe;
	const std::vector<double> values({ 0.0, 0.625 });

	CHECK_THROWS(TimeSerializerType(4, 1.25, 2, aoe, "Expression", values));
}
