#include <catch2/catch_test_macros.hpp>

#include <PsimagLite/InputPath.hpp>

#include <filesystem>
#include <fstream>

TEST_CASE("InputPath operations", "[InputPath]")
{
	PsimagLite::InputPath input_path;
	std::string           test_file = "test_input_path.txt";
	std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "dmrgpp_test_dir";
	std::filesystem::create_directories(temp_dir);

	std::filesystem::path file_path_in_temp = temp_dir / test_file;

	SECTION("Initial path is empty and finds current directory file")
	{
		std::ofstream ofs(test_file);
		ofs << "test";
		ofs.close();

		REQUIRE(std::filesystem::exists(test_file));
		auto found = input_path.findFirst(test_file);
		CHECK(found == std::filesystem::current_path() / test_file);

		std::filesystem::remove(test_file);
	}

	SECTION("Push new path and find file there")
	{
		std::ofstream ofs(file_path_in_temp);
		ofs << "test";
		ofs.close();

		input_path.push(temp_dir);
		auto found = input_path.findFirst(test_file);
		CHECK(found == file_path_in_temp.string());

		std::filesystem::remove(file_path_in_temp);
	}

	SECTION("Throws if file not found")
	{
		CHECK_THROWS_AS(input_path.findFirst("non_existent_file.txt"), std::runtime_error);
	}

	SECTION("LIFO behavior (last pushed path is searched first)")
	{
		std::filesystem::path temp_dir2
		    = std::filesystem::temp_directory_path() / "dmrgpp_test_dir2";
		std::filesystem::create_directories(temp_dir2);
		std::filesystem::path file_path_in_temp2 = temp_dir2 / test_file;

		std::ofstream ofs1(file_path_in_temp);
		ofs1 << "test1";
		ofs1.close();

		std::ofstream ofs2(file_path_in_temp2);
		ofs2 << "test2";
		ofs2.close();

		input_path.push(temp_dir);
		input_path.push(temp_dir2);

		auto found = input_path.findFirst(test_file);
		CHECK(found == file_path_in_temp2.string());

		std::filesystem::remove(file_path_in_temp);
		std::filesystem::remove(file_path_in_temp2);
		std::filesystem::remove(temp_dir2);
	}

	std::filesystem::remove_all(temp_dir);
}
